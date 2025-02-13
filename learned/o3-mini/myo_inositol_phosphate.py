"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
Definition: An inositol phosphate in which the inositol component has myo-configuration.
This improved heuristic determines if a given SMILES corresponds to a plain myo-inositol phosphate,
rejecting molecules where the inositol ring is acylated (i.e. having extra carbon chains outside
the six-membered ring).
"""

from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a plain myo-inositol phosphate based on its SMILES string.
    
    The procedure is:
      1. Parse the SMILES.
      2. Ensure at least one phosphorus atom is present.
      3. Identify a candidate six‐membered ring composed solely of carbons (the inositol core).
      4. Verify that each carbon in the ring has an oxygen substituent (not in the ring)
         and that none of these oxygens is “complex”: either they are hydroxyl groups or are 
         bound to a phosphorus atom. (If an oxygen is additionally bound to a carbon chain, we reject.)
      5. Confirm that every atom in the candidate ring has an explicit chiral tag.
      6. Crucially, ensure that no additional carbon atoms are present outside the inositol ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is classified as a myo-inositol phosphate.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Must contain a phosphorus atom.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found – cannot be an inositol phosphate"
    
    # Step 2: Identifying a candidate six-membered carbon ring (inositol core)
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    for ring in ring_info.AtomRings():
        # A candidate inositol ring is a 6-membered ring and every atom is carbon.
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            candidate_ring = ring
            break
    if candidate_ring is None:
        return False, "No six‐membered carbon ring found that could represent an inositol core"
    
    # Step 3: Reject if there are any extra carbon atoms outside the candidate ring.
    extra_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in candidate_ring]
    if extra_carbons:
        return False, "Found extra carbon atoms outside the inositol core – likely acylation present"
    
    # Step 4: Check stereochemistry on every atom in the candidate ring.
    from rdkit.Chem.rdchem import ChiralType
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
            return False, f"Inositol core atom at index {idx} lacks explicit stereochemistry"
    
    phosphate_found = False  # flag to ensure at least one substituent is a phosphate.
    
    # Step 5: For each carbon in the candidate ring, examine its substituents.
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        # Find neighbors that are not in the ring.
        oxygens = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in candidate_ring and nbr.GetAtomicNum() == 8]
        if not oxygens:
            return False, f"Inositol ring atom at index {idx} lacks any oxygen substituent"
        for oxy in oxygens:
            # Get substituents on oxygen excluding the inositol carbon.
            ext_neighbors = [nbr for nbr in oxy.GetNeighbors() if nbr.GetIdx() != idx]
            # Case 1: Oxygen bound to a phosphorus (phosphate substituent)
            if any(nbr.GetAtomicNum() == 15 for nbr in ext_neighbors):
                phosphate_found = True
                # Further check: ensure that the oxygen attached to phosphorus does not lead to extra carbons.
                for nbr in ext_neighbors:
                    if nbr.GetAtomicNum() == 15:
                        # Now check all neighbors of this phosphorus (except our oxygen)
                        for p_nbr in nbr.GetNeighbors():
                            if p_nbr.GetIdx() == oxy.GetIdx():
                                continue
                            # If any neighbor of phosphorus is carbon then it is not a plain phosphate.
                            if p_nbr.GetAtomicNum() == 6:
                                return False, (f"Phosphate group attached to inositol ring at atom index {idx} "
                                               "appears to be acylated (contains carbon bound to phosphorus)")
                continue  # oxygen is an acceptable phosphate substituent
            # Case 2: Oxygen not bound to phosphorus – should be a simple hydroxyl group.
            # We check that the oxygen is not attached to any carbon (which would indicate acylation).
            if any(nbr.GetAtomicNum() == 6 for nbr in ext_neighbors):
                return False, (f"Oxygen substituent on inositol ring at atom index {idx} is linked to a carbon chain; "
                               "suggesting acylation rather than a plain hydroxyl/phosphate")
            # If no objection, we accept this oxygen (hydroxyl).
    
    if not phosphate_found:
        return False, "Candidate inositol core does not have any phosphate substitution"
    
    return True, ("Molecule consists solely of a six-membered inositol core (all carbons stereochemically defined) "
                  "with only hydroxyl or simple phosphate substituents and no extra carbon chains – consistent with myo-inositol phosphate")
    
# Example test (using one positive from the prompt)
if __name__ == '__main__':
    test_smiles = "O[C@@H]1[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H]1O"  # 1D-myo-inositol 6-phosphate
    result, reason = is_myo_inositol_phosphate(test_smiles)
    print("Test molecule:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)
"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
Definition: An inositol phosphate in which the inositol component has myo-configuration.
This program heuristically determines if a given SMILES corresponds to a myo-inositol phosphate,
rejecting molecules where the inositol ring is acylated (as in phosphoinositides).
"""

from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
    The procedure is as follows:
      1. Parse the SMILES string into an RDKit molecule.
      2. Ensure the molecule contains at least one phosphorus atom.
      3. Identify a candidate six-membered ring composed solely of carbon atoms (the inositol core).
      4. Ensure every carbon of the candidate ring is substituted with at least one oxygen.
         Each oxygen substituent must be "simple": either a hydroxyl (only bound to the ring carbon)
         or directly bound to a phosphorus atom. If any substituent shows signs of complex acylation
         (i.e. is attached to a carbon chain), the molecule is rejected.
      5. Verify that every atom in the candidate ring has explicit stereochemical assignment,
         which is expected in the myo-configuration.
      6. Confirm that at least one of the oxygen substituents is linked to a phosphorus atom,
         to qualify the molecule as a phosphate.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the heuristic finds a myo-inositol phosphate, otherwise False.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule contains at least one phosphorus atom.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found – cannot be an inositol phosphate"
    
    # Identify candidate six-membered ring made exclusively of carbons.
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            candidate_ring = ring
            break
    if candidate_ring is None:
        return False, "No six‐membered carbon ring found that could represent an inositol core"
    
    # Check that every atom in the candidate inositol ring has an explicit chiral tag.
    from rdkit.Chem.rdchem import ChiralType
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
            return False, f"Atom at index {idx} in the candidate inositol core lacks explicit stereochemistry"
    
    phosphate_found = False
    # For each carbon in the candidate ring, check substituents.
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        # Look for at least one oxygen substituent on this ring atom (not part of the ring)
        oxygen_found = False
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in candidate_ring:
                continue  # we skip atoms in the ring
            if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                # Check the substituent: get its neighbors (excluding the inositol carbon)
                external_neighbors = [x for x in nbr.GetNeighbors() if x.GetIdx() != idx]
    
                # If oxygen is bound to a phosphorus, it qualifies as a phosphate substituent.
                if any(x.GetAtomicNum() == 15 for x in external_neighbors):
                    oxygen_found = True
                    phosphate_found = True
                    continue
                # Otherwise, if the oxygen is only attached to the inositol carbon (a hydroxyl group)
                # then that is acceptable if it is not linked to any carbon (which would indicate an acyl chain).
                if not external_neighbors:
                    oxygen_found = True
                    continue
                # If the oxygen is attached to another atom, check that it is not a carbon (indicative of acylation).
                if any(x.GetAtomicNum() == 6 for x in external_neighbors):
                    return False, (f"Oxygen substituent on inositol carbon at index {idx} "
                                   "is linked to a carbon chain; indicates acylation (phosphoinositide) rather than a plain inositol phosphate")
                # Otherwise, if only very simple atoms (like additional oxygens) are attached, we accept.
                oxygen_found = True
        if not oxygen_found:
            return False, f"Inositol ring atom at index {idx} lacks any oxygen substituent"
    
    if not phosphate_found:
        return False, "Candidate inositol core does not appear to have a phosphate substitution"
    
    return True, "Molecule contains a six-membered (inositol) ring with only hydroxyl or phosphate substituents and explicit stereochemistry, consistent with myo-inositol phosphate"


# Example test cases (the examples given in the user prompt; note: many are long, so here we show one positive)
if __name__ == '__main__':
    # Example positive: 1D-myo-inositol 6-phosphate
    test_smiles = "O[C@@H]1[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H]1O"
    result, reason = is_myo_inositol_phosphate(test_smiles)
    print("Test molecule:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)
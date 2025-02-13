"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
Definition: An inositol phosphate in which the inositol component has myo-configuration.
This improved heuristic:
  • parses the molecule,
  • identifies a single six‐membered carbon ring (the inositol core) with no extra carbon atoms,
  • verifies that each ring carbon has an oxygen substituent (which may be a simple hydroxyl or be phosphorylated),
  • “extracts” the inositol core (by excluding any extra atoms beyond the oxygen substituents of the ring),
  • computes a canonical InChIKey for that core and compares it to a reference D–myo–inositol,
  • and finally ensures that at least one oxygen substituent in the original molecule is bound to phosphorus.
If all conditions hold, the molecule is classified as a myo–inositol phosphate.
"""

from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a plain myo-inositol phosphate based on its SMILES string.

    The procedure is:
      1. Parse the SMILES.
      2. Ensure at least one phosphorus atom is present.
      3. Identify a candidate six‐membered ring composed solely of carbons (the inositol core).
      4. Reject if any extra carbon atoms are present (which could indicate acylation).
      5. For each ring carbon, verify that it has at least one oxygen substituent.
      6. Record if any of those oxygens is linked to a phosphorus.
      7. Extract the candidate “core” fragment (the ring carbons plus the directly attached oxygens)
         and for each oxygen, ignore any extra attachments (i.e. if phosphorylated, consider only the O).
      8. Compute a canonical identifier (via InChIKey) for this minimal inositol core and compare it with a
         reference D–myo–inositol core.
         
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if classified as a plain myo–inositol phosphate, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Must contain a phosphorus atom somewhere.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found – cannot be an inositol phosphate"

    ring_info = mol.GetRingInfo()
    candidate_ring = None
    # Look for a 6-membered ring where every atom is carbon.
    for ring in ring_info.AtomRings():
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            candidate_ring = ring
            break
    if candidate_ring is None:
        return False, "No six‐membered carbon ring found that could represent an inositol core"
    
    # Reject if any extra carbon atoms (outside the candidate ring) are present.
    extra_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in candidate_ring]
    if extra_carbons:
        return False, "Found extra carbon atoms outside the inositol core – likely acylation present"

    # For each ring carbon, ensure it has at least one oxygen substituent (outside the ring).
    # Also record if any oxygen is attached to phosphorus.
    phosphate_found = False
    candidate_fragment_indices = set(candidate_ring)  # will include ring C & the directly attached O's
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        # Get neighbors that are oxygen and not in the ring.
        oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in candidate_ring]
        if not oxygen_neighbors:
            return False, f"Inositol ring carbon at index {idx} lacks an oxygen substituent"
        # Add these oxygen atoms to our fragment.
        for oxy in oxygen_neighbors:
            candidate_fragment_indices.add(oxy.GetIdx())
            # Check if this oxygen is bound (in the full molecule) to phosphorus
            for nbr in oxy.GetNeighbors():
                if nbr.GetAtomicNum() == 15:
                    phosphate_found = True
                    break
    if not phosphate_found:
        return False, "Candidate inositol core does not have any phosphate substitution"

    # Now, we build a sub-molecule (fragment) containing:
    #   - the candidate ring carbons, and
    #   - their directly attached oxygens.
    # In this fragment, if an oxygen was originally phosphorylated,
    # the phosphorus (and its extra atoms) will be omitted.
    frag = Chem.PathToSubmol(mol, list(candidate_fragment_indices))
    # Important: assign stereochemistry on the fragment to preserve chiral tags.
    Chem.AssignStereochemistry(frag, cleanIt=True, force=True)
    
    # Now, obtain a canonical InChIKey for the fragment.
    try:
        frag_inchi = Chem.MolToInchiKey(frag)
    except Exception as e:
        return False, f"Error computing InChIKey for the candidate inositol core: {str(e)}"
    
    # Our reference for D-myo-inositol (all -OH groups at the positions) is defined by:
    ref_smiles = "O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    if ref_mol is None:
        return False, "Could not parse reference myo-inositol structure"
    Chem.AssignStereochemistry(ref_mol, cleanIt=True, force=True)
    try:
        ref_inchi = Chem.MolToInchiKey(ref_mol)
    except Exception as e:
        return False, f"Error computing InChIKey for reference myo-inositol: {str(e)}"
    
    # Compare the candidate inchi key (our extracted core) with that of myo-inositol.
    if frag_inchi != ref_inchi:
        return False, ("Extracted inositol core stereochemistry does not match D-myo-inositol. "
                       "Possibly a different inositol isomer (e.g. scyllo) or an acylated variant.")
    
    return True, ("Molecule consists solely of a six-membered inositol core (with all stereocenters defined) "
                  "whose minimal core (after replacing phosphorylated oxygens with hydroxyls) matches that of D-myo-inositol, "
                  "and at least one oxygen is phosphorylated – consistent with a myo-inositol phosphate.")

# Example test (using one positive from the prompt)
if __name__ == '__main__':
    # 1D-myo-inositol 6-phosphate example:
    test_smiles = "O[C@@H]1[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H]1O"
    result, reason = is_myo_inositol_phosphate(test_smiles)
    print("Test molecule:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)
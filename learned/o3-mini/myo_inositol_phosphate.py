"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
Definition: An inositol phosphate in which the inositol component has myo-configuration.
This program heuristically determines if a given SMILES corresponds to a myo-inositol phosphate.
"""

from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
    The procedure is as follows:
      1. Parse the SMILES string into an RDKit molecule.
      2. Ensure the molecule contains at least one phosphorus atom.
      3. Search for a candidate six-membered ring made solely of carbons (the inositol core).
      4. Verify that each carbon in the ring has at least one external oxygen substituent
         (corresponding to the hydroxyls or phosphate attached to the inositol ring).
      5. Confirm that at least one oxygen substituent from this ring is connected to a phosphorus atom.
      
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if the heuristic finds a myo-inositol phosphate, otherwise False.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule contains at least one phosphorus atom.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found – cannot be an inositol phosphate"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    # Look for a six-membered ring whose atoms are all carbon.
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                candidate_ring = ring
                break
    if candidate_ring is None:
        return False, "No six‐membered carbon ring found that could represent an inositol core"
    
    # For each carbon in the candidate ring check that it has at least one external oxygen substituent.
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        has_oxygen_substituent = False
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in candidate_ring and nbr.GetAtomicNum() == 8:
                has_oxygen_substituent = True
                break
        if not has_oxygen_substituent:
            return False, f"Atom at index {idx} in the candidate inositol core lacks an oxygen substituent"
    
    # Check that at least one of the oxygen substituents is linked to a phosphorus atom.
    phosphate_found = False
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in candidate_ring and nbr.GetAtomicNum() == 8:
                # Check neighbors of oxygen for phosphorus.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 15:
                        phosphate_found = True
                        break
                if phosphate_found:
                    break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "Candidate inositol core does not appear to have a phosphate substitution"

    # If all tests are passed, we conclude the molecule is a myo-inositol phosphate.
    return True, "Molecule contains a six-membered (inositol) ring with oxygen substituents and at least one phosphate group, consistent with myo-inositol phosphate"
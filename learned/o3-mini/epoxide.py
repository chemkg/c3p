"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide – Any cyclic ether in which the oxygen atom forms part of a 3‐membered ring.
For our purposes we only classify an epoxide if the 3‐membered oxirane ring is isolated 
(i.e. none of the three ring atoms belong to any other ring).
"""

from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide (isolated oxirane) based on its SMILES string.
    Our detection requires a 3-membered ring containing two tetrahedral carbons (sp³) 
    and one oxygen (with exactly 2 heavy-atom bonds) where none of the ring atoms participate 
    in any additional ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an epoxide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Precompute ring membership counts for isolation check
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if not atom_rings:
        return False, "No rings found in molecule"
    
    # Build a dictionary: for each atom index, count in how many rings it is found
    atom_ring_count = {atom.GetIdx(): 0 for atom in mol.GetAtoms()}
    for ring in atom_rings:
        for idx in ring:
            atom_ring_count[idx] += 1

    # Define a SMARTS pattern for an oxirane group:
    # – two tetrahedral carbons (C;X4) in a ring of size 3
    # – one oxygen with two heavy neighbors (O;X2) in the same 3‐membered ring.
    # The ;R3 modifier ensures that the atom is in a ring of exactly three atoms.
    epoxide_smarts = "[C;X4;R3][O;X2;R3][C;X4;R3]"
    query = Chem.MolFromSmarts(epoxide_smarts)
    if query is None:
        return False, "Error in SMARTS pattern"

    # Find all matches to the SMARTS 
    matches = mol.GetSubstructMatches(query)
    for match in matches:
        # match is a tuple of three atom indices
        # confirm that these atoms form one of the 3-membered rings in the molecule.
        # use the isolation criterion: each atom must appear in exactly one ring.
        if all(atom_ring_count[idx] == 1 for idx in match):
            # Additionally, check that the oxygen (middle atom in our pattern) has exactly 2 neighbors.
            # (This ensures it is not further substituted in a way that might alter its reactivity.)
            oxygen_atom = mol.GetAtomWithIdx(match[1])
            # Use GetDegree (ignoring implicit hydrogens) 
            if oxygen_atom.GetDegree() != 2:
                continue
            return True, "Found an isolated 3-membered cyclic ether (epoxide) group"
    
    return False, "Does not contain an isolated 3-membered cyclic ether (epoxide) group"


# Example usage (this block may be removed or commented out as needed)
if __name__ == "__main__":
    test_examples = [
        ("ClCC1CO1", "epichlorohydrin (should be True)"),
        ("C1CC1", "Cyclopropane (should be False)"),
        ("O1C(C)C1", "Epoxide-like ring (should be True if isolated)")
    ]
    for smi, desc in test_examples:
        result, reason = is_epoxide(smi)
        print(f"SMILES: {smi} | {desc} --> {result}, {reason}")
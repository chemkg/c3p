"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: Glycerophosphocholine
Definition: The glycerol phosphate ester of a phosphocholine. A nutrient with many different roles in human health.
A key structural feature is the phosphocholine unit – a phosphate group (P) attached to an oxygen which in turn is linked to a choline moiety (OCC[N+](C)(C)C).
This program uses RDKit to check for that connectivity.
"""
from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is classified as glycerophosphocholine based on its SMILES string.
    The algorithm requires that the molecule contains a phosphocholine moiety, 
    i.e. an oxygen-linked choline group (OCC[N+](C)(C)C) attached to a phosphorus atom.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a glycerophosphocholine, False otherwise
        str: Explanation for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the choline fragment: –OCC[N+](C)(C)C
    choline_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")
    if choline_pattern is None:
        return False, "Error creating choline pattern"
    
    # Search for choline fragment(s) in molecule
    choline_matches = mol.GetSubstructMatches(choline_pattern)
    if not choline_matches:
        return False, "Choline fragment not found"
    
    # For each match, check if the oxygen (first atom in the pattern) is attached to a phosphorus atom.
    for match in choline_matches:
        # match[0] corresponds to the oxygen in the choline moiety.
        oxygen_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:  # Phosphorus has atomic number 15
                # Additional check: Does the P atom have at least one double-bonded oxygen? (Optional check.)
                p_atom = neighbor
                has_double_bond_O = any(nbr for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and
                                        mol.GetBondBetweenAtoms(p_atom.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble() == 1)
                if has_double_bond_O:
                    return True, "Found phosphocholine moiety: choline group connected via oxygen to a phosphate group."
                else:
                    # Even if the P does not have a confirmed double-bonded oxygen,
                    # the connectivity can be indicative.
                    return True, "Found phosphocholine moiety: choline group connected via oxygen to phosphorus."
    
    return False, "No phosphorus atom found attached to the choline group"
    
# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES for glycerophosphocholine-like structure
    test_smiles = "P(OC[C@@H](COC(CCCCCCCCCCCCCCCCC)=O)OC(=O)CCCC)(=O)(OCC[N+](C)(C)C)[O-]"
    result, reason = is_glycerophosphocholine(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
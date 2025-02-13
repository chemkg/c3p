"""
Classifies: CHEBI:24128 furanocoumarin
"""
#!/usr/bin/env python
"""
Classifies: Any furochromene that consists of a furan ring fused with a coumarin
(i.e. furanocoumarins)
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a coumarin (benzopyrone) with a fused furan ring.
    The fusion is identified by searching for a coumarin core (lactone-fused benzene)
    and an aromatic furan ring that share at least two atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a furanocoumarin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the coumarin core.
    # This pattern represents benzopyrone: a benzene ring fused with an α-pyrone.
    # It is defined as “c1ccc2C(=O)Oc2c1” where:
    #   - The lactone is formed by C(=O)O
    #   - The rings share two atoms (positions 1 and 2)
    coumarin_smarts = "c1ccc2C(=O)Oc2c1"
    coumarin_pattern = Chem.MolFromSmarts(coumarin_smarts)
    if coumarin_pattern is None:
        return False, "Invalid coumarin SMARTS pattern"
    
    # Define a SMARTS pattern for an aromatic furan ring.
    # Furan is a five-membered ring with one oxygen.
    furan_smarts = "c1ccoc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS pattern"
    
    # Find all substructure matches for the coumarin core and the furan ring.
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    if not coumarin_matches:
        return False, "Coumarin (benzopyrone) core not found"
    
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "Furan ring not found"
    
    # Check if any coumarin match and furan match are fused, e.g.
    # if they share at least two atoms. This is a heuristic check for fusion.
    fusion_found = False
    for c_match in coumarin_matches:
        set_c = set(c_match)
        for f_match in furan_matches:
            set_f = set(f_match)
            # Intersection: the atoms that are in both patterns
            common_atoms = set_c.intersection(set_f)
            if len(common_atoms) >= 2:
                fusion_found = True
                break
        if fusion_found:
            break
    
    if not fusion_found:
        return False, "No fused furan and coumarin (benzopyrone) substructures found"
    
    # If we've reached here, we have found a coumarin core fused with a furan ring.
    return True, "Contains a coumarin core fused with a furan ring (furanocoumarin)"
    
# Example usage:
if __name__ == "__main__":
    test_smiles = "COc1cc2occc2c2oc(=O)ccc12"  # Example: Isobergapten
    result, reason = is_furanocoumarin(test_smiles)
    print(result, ":", reason)
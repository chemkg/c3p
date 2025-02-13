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
    This function searches for the coumarin core (using an updated SMARTS pattern to capture
    aromatic coumarin cores) and an aromatic furan ring, then checks whether they are fused 
    (sharing at least two atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a furanocoumarin, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the coumarin core.
    # Coumarin (benzopyrone) is now represented with an aromatic pattern:
    # "c1ccc2oc(=O)cc2c1" where:
    #    - The aromatic benzene ring (c1ccc...c1)
    #    - The fused lactone ring (oc(=O)cc) where the carbonyl is included.
    coumarin_smarts = "c1ccc2oc(=O)cc2c1"
    coumarin_pattern = Chem.MolFromSmarts(coumarin_smarts)
    if coumarin_pattern is None:
        return False, "Invalid coumarin SMARTS pattern"
    
    # Define a SMARTS pattern for an aromatic furan ring.
    # Furan is a five-membered aromatic ring containing one oxygen.
    furan_smarts = "c1ccoc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS pattern"
    
    # Find all matches for the coumarin core
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    if not coumarin_matches:
        return False, "Coumarin (benzopyrone) core not found"
    
    # Find all matches for the furan ring
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "Furan ring not found"
    
    # Check whether any coumarin match and furan match are fused
    # Fusion is assumed if at least two atoms are shared between the matches.
    fusion_found = False
    for c_match in coumarin_matches:
        set_c = set(c_match)
        for f_match in furan_matches:
            set_f = set(f_match)
            common_atoms = set_c.intersection(set_f)
            if len(common_atoms) >= 2:
                fusion_found = True
                break
        if fusion_found:
            break
    
    if not fusion_found:
        return False, "No fused furan and coumarin (benzopyrone) substructures found"
    
    # If both substructures are found and are fused, then classify as a furanocoumarin.
    return True, "Contains a coumarin core fused with a furan ring (furanocoumarin)"

# Example usage:
if __name__ == "__main__":
    # Example: Isobergapten SMILES: COc1cc2occc2c2oc(=O)ccc12
    test_smiles = "COc1cc2occc2c2oc(=O)ccc12"
    result, reason = is_furanocoumarin(test_smiles)
    print(result, ":", reason)
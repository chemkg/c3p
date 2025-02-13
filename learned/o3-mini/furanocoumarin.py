"""
Classifies: CHEBI:24128 furanocoumarin
"""
#!/usr/bin/env python
"""
Classifies: Any furochromene that consists of a furan ring fused with a coumarin (furanocoumarin)
The function is_furanocoumarin uses improved SMARTS patterns for the coumarin (benzopyrone) core
and verifies that an aromatic furan ring is fused (sharing >=2 atoms) with the coumarin core.
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a coumarin (benzopyrone) with a fused furan ring.
    This function searches for the coumarin core (using one of two SMARTS variations to be more tolerant)
    and for an aromatic furan ring. Then it checks whether a furan match shares at least two atoms 
    with any coumarin match, which is our proxy for fusion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a furanocoumarin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define two SMARTS patterns for the coumarin (benzopyrone) core.
    # The first is the classic pattern and the second a slight variation to catch substituted or fused cases.
    coumarin_smarts_list = [
        "c1ccc2oc(=O)cc2c1",   # common representation
        "c1ccc2c(c1)oc(=O)cc2"  # variation (ring order swapped)
    ]
    
    coumarin_matches = []
    for smarts in coumarin_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            coumarin_matches.extend(matches)
    if not coumarin_matches:
        return False, "Coumarin (benzopyrone) core not found"
    
    # Define a SMARTS pattern for an aromatic furan ring.
    # Furan is a five-membered ring with one oxygen atom.
    furan_smarts = "c1ccoc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS pattern"
    
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "Furan ring not found"
    
    # Check for fusion: require that at least one furan match and one coumarin match share
    # two or more atoms. This is our proxy for the rings being fused.
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
        return False, "Coumarin core and furan ring not properly fused"
    
    return True, "Contains a coumarin (benzopyrone) core fused with a furan ring (furanocoumarin)"


# Example usage:
if __name__ == "__main__":
    # Test with one of the known furanocoumarins - for instance, Isobergapten
    test_smiles = "COc1cc2occc2c2oc(=O)ccc12"
    result, reason = is_furanocoumarin(test_smiles)
    print(result, ":", reason)
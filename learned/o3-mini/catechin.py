"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechin class – Members of the class of hydroxyflavan that have a flavan-3-ol skeleton and its substituted derivatives.
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin class based on its SMILES string.
    The catechin (flavan-3-ol) core is roughly represented as the 2-phenyl-3,4-dihydro-2H-chromene-3-ol scaffold.
    This minimal pattern consists of an aromatic ring (A-ring) fused to an oxygen-containing heterocycle (C-ring)
    that carries at least one hydroxyl (at C3) and is substituted by an extra aromatic ring (B-ring).
    
    Note: Due to the many substitutions found in natural catechins and their derivatives, 
    this SMARTS pattern is a simplified heuristic and may not match every derivative.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a catechin (flavan-3-ol), False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define a SMARTS pattern that roughly captures the flavan-3-ol core.
    # This pattern “c1ccc2OC(C(O)C2)c1” looks for:
    #   - an aromatic ring (c1ccc...c1)
    #   - fused to an oxygen-containing heterocycle (OC(C(O)C2)) 
    #   - with an extra aromatic group attached (the C-ring bears a substituted carbon from which the second aryl is connected).
    # (In practice, this simplified pattern is a heuristic.)
    core_smarts = "c1ccc2OC(C(O)C2)c1"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return None, None  # Should not occur
    
    # Check for the substructure match of the core scaffold.
    if mol.HasSubstructMatch(core_pattern):
        return True, "Molecule contains a core flavan-3-ol scaffold characteristic of catechins."
    else:
        return False, "Molecule does not appear to contain the required flavan-3-ol (catechin) scaffold."

# Example usage:
if __name__ == "__main__":
    # Sample SMILES for (–)-catechin (one of the examples)
    test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"
    is_cat, reason = is_catechin(test_smiles)
    print("Is catechin?", is_cat)
    print("Reason:", reason)
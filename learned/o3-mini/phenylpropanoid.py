"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: phenylpropanoid compounds
Definition: Any organic aromatic compound with a structure based on a
phenylpropane (C6–C3) skeleton. This class includes naturally occurring
phenylpropanoid esters, flavonoids, anthocyanins, coumarins and many small
phenolic molecules as well as their semi-synthetic and synthetic analogues.
Phenylpropanoids are also precursors of lignin.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    
    The heuristic checks:
      - That a benzene ring is present.
      - That it also contains a C6–C3 fragment via one of several SMARTS patterns:
          · A saturated three‐carbon chain attached to a benzene ring.
          · An unsaturated three‐carbon chain (cinnamic type) attached to benzene.
          · A cinnamic acid fragment.
          · A coumarin substructure.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a phenylpropanoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize molecule (computes aromaticity etc.)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Basic check: must contain at least one benzene ring.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene ring detected, hence not a phenylpropanoid"

    # Define a list of (description, SMARTS) that represent common phenylpropanoid fragments.
    patterns = [
        ("Benzene with a saturated propyl chain (C6–C3)", "c1ccccc1CCC"),
        ("Benzene with an unsaturated (allylic) chain (C6–C3)", "c1ccccc1C=CC"),
        ("Cinnamic acid type fragment", "c1ccccc1C=CC(=O)[O-]"),  # note: acid may be deprotonated too
        ("Coumarin structure", "O=C1Oc2ccccc2C1")
    ]
    
    # Try each pattern and if any matches, we classify as phenylpropanoid.
    for desc, smarts in patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip if SMARTS did not compile
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches pattern: {desc}"
    
    # Optionally, check for a flavonoid-like core (C6–C3–C6) as an additional heuristic.
    # A simple flavonoid core: two aromatic rings linked via a heterocycle.
    flavonoid_pattern = Chem.MolFromSmarts("c1ccc2c(c1)OC(=O)[C@H]1c(O)ccc(O)c1C2") 
    if flavonoid_pattern and mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Matches a flavonoid-type substructure (C6–C3–C6)"
    
    return False, "No phenylpropanoid-characteristic substructure found"

# Example usage:
if __name__ == "__main__":
    # A few test cases (the SMILES given in the prompt are complex, so here is a simpler one)
    test_smiles = [
        "c1ccccc1CCC",  # benzene with a propyl chain (should be phenylpropanoid-like)
        "CCO",         # ethanol (should not match)
    ]
    
    for sm in test_smiles:
        result, reason = is_phenylpropanoid(sm)
        print(f"SMILES: {sm}\nClassified as phenylpropanoid? {result}\nReason: {reason}\n")
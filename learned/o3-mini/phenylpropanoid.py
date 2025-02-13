"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: phenylpropanoid compounds
Definition: Any organic aromatic compound with a structure based on a 
phenylpropane (C6–C3) skeleton. This class includes naturally occurring 
phenylpropanoid esters, flavonoids, anthocyanins, coumarins and many small 
phenolic molecules as well as their semi‐synthetic and synthetic analogues.
Phenylpropanoids are also precursors of lignin.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    
    The heuristic implemented here requires:
      - That the molecule contains a benzene ring.
      - That one of several phenylpropanoid-characteristic substructures is present.
        These include:
          · Coumarin substructure
          · A cinnamic acid–like fragment (benzene with a C=CC(=O)[O] fragment)
          · Benzene attached to an unsaturated 3-carbon chain (cinnamyl type)
          · Benzene attached to a saturated propyl chain
          · A flavanone core
          
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a phenylpropanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule (to compute aromaticity, assign bond orders, etc.)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Basic filter: Check that at least one benzene ring is present
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene ring detected, so not a phenylpropanoid"
    
    # List of SMARTS patterns along with human‐readable descriptions.
    # These patterns have been refined to capture the C6–C3 scaffold in various contexts.
    patterns = [
        ("Coumarin structure", "O=C1Oc2ccccc2C1"),
        # The cinnamic acid fragment: benzene attached to an unsaturated three-carbon chain bearing a carbonyl and an oxygen.
        ("Cinnamic acid fragment", "c1ccccc1C=CC(=O)[O;H0,-1]"),
        ("Benzene with unsaturated (cinnamyl) chain", "c1ccccc1C=CC"),
        ("Benzene with saturated (propyl) chain", "c1ccccc1CCC"),
        # A simple flavanone core represented by two aromatic rings linked via a heterocycle.
        ("Flavanone core", "c1cc(c(cc1)O)C2=CC(=O)OC3=CC=CC=C23")
    ]
    
    # Check each pattern: if any match, then classify as a phenylpropanoid.
    for desc, smarts in patterns:
        frag = Chem.MolFromSmarts(smarts)
        if frag is None:
            continue  # If a SMARTS pattern fails to compile, skip it.
        if mol.HasSubstructMatch(frag):
            return True, f"Matches pattern: {desc}"
    
    return False, "No phenylpropanoid-characteristic substructure found"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "c1ccccc1CCC",  # benzene with a saturated propyl chain (expected phenylpropanoid-like)
        "c1ccccc1C=CC", # benzene with an unsaturated chain (expected phenylpropanoid-like)
        "OC(=Cc1ccccc1)C([O-])=O",  # enol-phenylpyruvate (should not be classified)
        "CCO"  # ethanol (should not be classified)
    ]
    
    for sm in test_smiles:
        result, reason = is_phenylpropanoid(sm)
        print(f"SMILES: {sm}\nClassified as phenylpropanoid? {result}\nReason: {reason}\n")
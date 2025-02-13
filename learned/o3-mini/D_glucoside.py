"""
Classifies: CHEBI:35436 D-glucoside
"""
#!/usr/bin/env python3
"""
Classifies: D-glucoside (any glucoside in which the glycoside group is derived from D-glucose)
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    
    A D-glucoside is defined here as a molecule containing a glycosidic linkage derived from D-glucose.
    In many representations the sugar unit appears as a glucopyranose ring (six-membered ring with 5 carbons
    and 1 oxygen) substituted at the anomeric carbon via an O-glycosidic bond. Common depictions for the sugar
    moiety are (for example) one of the two patterns below (one for the alpha-anomer and one for the beta-anomer).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a D-glucoside, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for an O-glycosidically linked D-glucopyranosyl unit.
    # These patterns try to capture a common depiction of the sugar unit.
    #
    # The sugar ring is a six-membered ring with one ring oxygen and five carbons.
    # In a glucoside the anomeric hydroxyl of glucose is replaced by an oxygen-linked substituent.
    # The patterns below start with an extra O (the bridging, glycosidic oxygen) attached to a chiral carbon
    # in a ring that has a CH2OH substituent (denoted by CO) at one position.
    #
    # Note: Because of variability in depictions, these SMARTS may miss some true positives or include some false ones.
    alpha_pattern = Chem.MolFromSmarts("O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    beta_pattern  = Chem.MolFromSmarts("O[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O")
    
    # Search for either pattern in the molecule.
    if mol.HasSubstructMatch(alpha_pattern):
        return True, "Found an alpha-D-glucopyranosyl moiety (D-glucoside)"
    if mol.HasSubstructMatch(beta_pattern):
        return True, "Found a beta-D-glucopyranosyl moiety (D-glucoside)"
    
    return False, "No D-glucoside (D-glucopyranosyl) moiety found in the molecule"

# For testing purposes:
if __name__ == '__main__':
    # Example SMILES strings (taken from the provided examples):
    test_smiles = [
        # beta-D-glucosyl-(1<->1')-N-[(17Z)-hexacosenoyl]sphinganine
        "CCCCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCCCCCCCCC\\C=C/CCCCCCCC",
        # apigenin 7-O-beta-D-glucoside
        "C1=2C(C(C=C(O1)C=3C=CC(O)=CC3)=O)=C(C=C(C2)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O",
    ]
    
    for s in test_smiles:
        result, reason = is_D_glucoside(s)
        print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")
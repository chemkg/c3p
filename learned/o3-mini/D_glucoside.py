"""
Classifies: CHEBI:35436 D-glucoside
"""
#!/usr/bin/env python3
"""
Classifies: D-glucoside (any glucoside in which the glycoside group is derived from D-glucose)

This program attempts to detect a D-glucopyranosyl moiety by searching for a six-membered ring (pyranose)
with one ring oxygen, a CH2OH substituent at the C6 position, and an extra oxygen at the anomeric carbon.
Two SMARTS patterns are provided to cover both beta and alpha anomers.
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.

    A D-glucoside is defined as a molecule containing a glucosyl unit derived from D-glucose,
    typically depicted as a pyranose ring (six-membered with five carbons and one oxygen) that includes
    a CH2OH substituent at the C6 position and a glycosidic linkage via the anomeric oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a D-glucoside, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # -------------------------------
    # The idea is to capture a D-glucopyranosyl moiety.
    # We define two SMARTS patterns which include:
    # 1. [O;D2]           : a bridging glycosidic oxygen (exactly two connections)
    # 2. [C@H] or [C@@H]1 : the chiral anomeric carbon starting a ring labelled '1'
    # 3. A pyran ring completion with the presence of appropriate hydroxyl groups,
    #    and a terminal CH2OH (CO) group corresponding to C6.
    #
    # The two patterns below account for beta and alpha configurations at the anomeric carbon.
    # -------------------------------
    beta_smarts = "[O;D2][C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1CO"
    beta_pattern = Chem.MolFromSmarts(beta_smarts)
    
    alpha_smarts = "[O;D2][C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1CO"
    alpha_pattern = Chem.MolFromSmarts(alpha_smarts)
    
    # Search for the beta pattern first, then the alpha pattern.
    if mol.HasSubstructMatch(beta_pattern):
        return True, "Found a beta-D-glucopyranosyl moiety (D-glucoside)"
    if mol.HasSubstructMatch(alpha_pattern):
        return True, "Found an alpha-D-glucopyranosyl moiety (D-glucoside)"
    
    return False, "No D-glucoside (D-glucopyranosyl) moiety found in the molecule"


# The following block can be used for testing purposes.
if __name__ == '__main__':
    # Examples of SMILES strings known to be D-glucosides.
    test_smiles = [
        # beta-D-glucosyl-(1<->1')-N-[(17Z)-hexacosenoyl]sphinganine:
        "CCCCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCCCCCCCCC\\C=C/CCCCCCCC",
        # apigenin 7-O-beta-D-glucoside:
        "C1=2C(C(C=C(O1)C=3C=CC(O)=CC3)=O)=C(C=C(C2)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O",
    ]
    
    for s in test_smiles:
        result, reason = is_D_glucoside(s)
        print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")
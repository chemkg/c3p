"""
Classifies: CHEBI:35436 D-glucoside
"""
Below is one possible improved solution. In our reasoning we observed that the previous attempt simply looked for two very narrowly defined SMARTS strings that were too strict (and in some cases “over‐matching” fragments already present in complex multimers) while missing many representations. In our new approach we try to capture a more “generic” D‐glucopyranosyl fragment. For example, we require (a) that the sugar ring is six–membered and contains one ring oxygen, (b) that one of the ring carbons bears the CH2OH substituent (characteristic of a D‐glucose moiety), and (c) that an additional oxygen is bound at the anomeric centre (the glycosidic linkage) rather than a hydroxyl group. (Note that many sugar depictions vary, so no single SMARTS will be 100% perfect.) 

In the code below we define two SMARTS patterns – one for the beta– and one for the alpha–anomer – that demand that an extra (bridging) O is attached to the anomeric carbon of a glucopyranosyl ring. The SMARTS uses explicit chirality tags and also checks that the substituent at the “6–position” is a CH2OH moiety. Finally, we search the molecule for either pattern and return a positive classification only if one is found. (If the sugar is depicted in a nonstandard way or has extra modifications, it may still be missed.) 

Below is the complete Python program:

------------------------------------------------

#!/usr/bin/env python3
"""
Classifies: D-glucoside (any glucoside in which the glycoside group is derived from D-glucose)

This improved program attempts to detect a D-glucopyranosyl moiety by searching for 
a six-membered ring (i.e. pyranose) with one ring oxygen, a CH2OH substituent at the position
that (in D-glucose) corresponds to C6, and an extra oxygen bound at the anomeric carbon 
(the glycosidic oxygen). Two SMARTS patterns (for beta and alpha anomers) are provided.
Due to variability in depictions, the method may still yield false positives/negatives.
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.

    A D-glucoside is defined here as a molecule containing a glucosyl unit derived from D-glucose,
    typically found as a pyranose ring (six-membered ring with five carbons and one oxygen) having a CH2OH
    substituent at the C6 position, with the anomeric carbon (C1) linked via an oxygen to a substituent 
    (the glycoside linkage). Two SMARTS patterns (for alpha and beta stereochemistry) are defined.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a D-glucoside, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ------------------------------
    # Improved SMARTS explanations:
    #
    # 1. We require an extra oxygen (the glycosidic oxygen) attached to the anomeric carbon.
    #    We write this as [O;D2] meaning an oxygen with 2 bonds.
    #
    # 2. The ring should be a pyranose (six-membered) ring; we require five chiral atoms,
    #    four of which are substituted with –OH, and one (the 6–position) bearing a CH2OH group.
    #
    # 3. The two patterns below roughly capture the beta and alpha stereochemistry at the anomeric center.
    #    (Exact chirality can be notationally inconsistent between sources so these are best–guesses.)
    #
    # The patterns below are not perfect. For example, they assume that the –CH2OH group attached to the 
    # appropriate ring carbon is depicted explicitly and not further substituted.
    # ------------------------------

    # SMARTS for a beta-D-glucopyranosyl moiety:
    # The fragment defined is:
    #   [O;D2]           : a bridging glycosidic oxygen (2 connections)
    #   [C@H]1           : anomeric carbon with specified chirality, beginning ring 1
    #   O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1CO : rest of the six-membered ring with a CH2OH substituent (CO) at the end.
    beta_smarts = "[O;D2][C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1CO"
    beta_pattern = Chem.MolFromSmarts(beta_smarts)
    
    # SMARTS for an alpha-D-glucopyranosyl moiety:
    # Similar to above but with the anomeric carbon chirality inverted.
    alpha_smarts = "[O;D2][C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1CO"
    alpha_pattern = Chem.MolFromSmarts(alpha_smarts)
    
    # Now search for either pattern in the input molecule.
    # The match (if any) should correspond to a glycosidically linked D-glucose unit.
    if mol.HasSubstructMatch(beta_pattern):
        return True, "Found a beta-D-glucopyranosyl moiety (D-glucoside)"
    if mol.HasSubstructMatch(alpha_pattern):
        return True, "Found an alpha-D-glucopyranosyl moiety (D-glucoside)"
    
    return False, "No D-glucoside (D-glucopyranosyl) moiety found in the molecule"

# For testing purposes:
if __name__ == '__main__':
    # Some test SMILES strings (excerpts from provided examples):
    test_smiles = [
        # beta-D-glucosyl-(1<->1')-N-[(17Z)-hexacosenoyl]sphinganine
        "CCCCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCCCCCCCCC\\C=C/CCCCCCCC",
        # apigenin 7-O-beta-D-glucoside
        "C1=2C(C(C=C(O1)C=3C=CC(O)=CC3)=O)=C(C=C(C2)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O",
    ]
    
    for s in test_smiles:
        result, reason = is_D_glucoside(s)
        print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")

------------------------------------------------

Notes:
• This approach makes several assumptions about sugar ring depiction and chirality.
• Despite these improvements, variability in chemical depictions means that some true D‐glucosides
  may be missed (false negatives) and some other fragments might spuriously match (false positives).
• Further improvement could include additional pattern variants or full sugar‐recognition algorithms.

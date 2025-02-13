"""
Classifies: CHEBI:35436 D-glucoside
"""
#!/usr/bin/env python3
"""
Classifies: D-glucoside 
Definition: A D-glucoside is any glucoside in which the glycoside group is derived from D-glucose.
This program attempts to detect a D-glucoside unit by first matching a relaxed SMARTS pattern
representing the connectivity of a glucopyranosyl moiety – namely a six‐membered ring containing one oxygen
(with an exocyclic CH2OH substituent) – and then verifying that the extracted fragment has the formula C6H11O5.
This fragment formula corresponds to a glucosyl unit (derived from D-glucose) after formation of a glycosidic bond.
"""

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    
    The detection is done in two steps:
      1. A relaxed SMARTS is used to find a pyranose‐like ring with connectivity expected
         for a glucosyl residue (one ring oxygen and a CH2OH group attached at the ring).
      2. For each match, the substructure is extracted and its formula is computed.
         For a glucoside derived from D-glucose, the glycosyl (or glucosyl) unit should have the formula C6H11O5.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a D-glucoside moiety is found, False otherwise.
        str: A message describing the outcome.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a relaxed SMARTS pattern to capture a glucosyl unit.
    # The pattern looks for:
    #   - An oxygen (the glycosidic oxygen) connected to
    #   - A ring (labelled with "1") with the pattern: C1 O C(O) C(O) C(O) C1 CO
    #
    # This roughly represents a glucopyranosyl moiety after glycosidic bond formation.
    # (The stereochemistry markers are omitted so as to capture more examples.)
    relaxed_smarts = "O[C]1OC(O)C(O)C(O)C1CO"
    pattern = Chem.MolFromSmarts(relaxed_smarts)
    if pattern is None:
        return False, "Failed to compile SMARTS pattern"
    
    # Look for substructure matches. We set useChirality=False to avoid missing matches
    # that do not have explicit chirality flags.
    matches = mol.GetSubstructMatches(pattern, useChirality=False)
    if not matches:
        return False, "No sugar moiety matching the expected glucosyl connectivity found"
    
    # For each matching fragment, extract it and check its calculated molecular formula.
    # A D-glucoside (after glycosidic bond formation) should have the formula C6H11O5.
    for match in matches:
        # Create a sub-molecule from the atoms in the match.
        submol = Chem.PathToSubmol(mol, list(match))
        formula = CalcMolFormula(submol)
        # Uncomment the following line to see details during debugging:
        # print(f"DEBUG: Found fragment formula: {formula}")
        if formula == "C6H11O5":
            return True, "Found a D-glucoside moiety (glucopyranosyl fragment with formula C6H11O5)"
    return False, "Found sugar-like fragment(s) but none matching the D-glucoside formula C6H11O5"

# For testing purposes (run only when invoked as the main program)
if __name__ == '__main__':
    # Here are some example SMILES (some from the user’s examples).
    test_smiles = [
        # Expected to be a D-glucoside
        "CCCCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCCCCCCCCC\\C=C/CCCCCCCC",
        "C1=2C(C(C=C(O1)C=3C=CC(O)=CC3)=O)=C(C=C(C2)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O",
        # Expected not to be a D-glucoside
        "ClC1=C(OC2OC(O)C(O)C(C2C(=O)O)O)C=CC(=C1)C3NC(=O)[C@H](N(C(=O)[C@@H](NC(=O)[C@@H](C)CC(=C[C@@H]([C@H](OC([C@H]3OC)=O)C)C)C)C)C)CC4=C(Cl)NC5=C4C=CC=C5"
    ]
    
    for s in test_smiles:
        result, reason = is_D_glucoside(s)
        print(f"SMILES: {s}\nClassification: {result}\nReason: {reason}\n")
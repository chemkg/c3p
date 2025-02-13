"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
Definition: Any D-galactoside having beta-configuration at its anomeric centre.
This program uses a heuristic SMARTS approach to look for a pyranose ring
with a glycosidic substituent at the anomeric carbon and with the expected stereochemistry.
Note: Carbohydrate isomers are difficult to identify computationally; if no match is found,
the function returns False with an appropriate reason.
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    Heuristics:
      - It must contain a galactopyranoside (six-membered sugar ring with one oxygen).
      - The anomeric carbon (the carbon attached via an -O- bond to an external residue)
        should have the appropriate beta stereochemical annotation.
      - Allow for modifications (e.g. acetylation) at positions on the ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a beta-D-galactoside, else False.
        str: A reason for the classification.
    """
    # Parse SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # We now define two SMARTS patterns that capture the beta-D-galactoside core.
    # The patterns below look for a six-membered ring (pyranose) that:
    #   - Contains 5 carbons and 1 oxygen.
    #   - Has an exocyclic substituent attached via an oxygen at the anomeric position.
    # The wildcard ([*:]) allows for many families of substituents.
    #
    # Note: The chiral tags ([C@@H] and [C@H]) are used to capture the relative stereochemistry.
    # Given differences in drawing/convention the beta annotation may be encoded either way,
    # so we offer two SMARTS variants.
    smarts1 = "[*:1]-O[C@@H]1O[C@H]([*:2])[C@H](O)[C@H](O)[C@H]1O"
    smarts2 = "[*:1]-O[C@H]1O[C@@H]([*:2])[C@@H](O)[C@H](O)[C@@H]1O"
    
    patt1 = Chem.MolFromSmarts(smarts1)
    patt2 = Chem.MolFromSmarts(smarts2)
    
    if patt1 is None or patt2 is None:
        # in case SMARTS could not be parsed (should not happen here)
        return None, None

    # Check if at least one of the SMARTS patterns is found as a substructure match.
    if mol.HasSubstructMatch(patt1) or mol.HasSubstructMatch(patt2):
        return True, "Detected beta-D-galactoside moiety (pyranose ring with beta glycosidic linkage)"
    else:
        return False, "No beta-D-galactoside moiety detected"

# Example test cases (uncomment to run)
# test_smiles = [
#     "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",  # methyl beta-D-galactoside
#     "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@H]1O"  # cyanidin 3-O-beta-D-galactoside
# ]
# for sm in test_smiles:
#     flag, reason = is_beta_D_galactoside(sm)
#     print(f"SMILES: {sm}\nMatch: {flag}\nReason: {reason}\n")
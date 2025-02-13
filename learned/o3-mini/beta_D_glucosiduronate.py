"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group of any beta-D-glucuronic acid.
Major species at pH 7.3.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    
    Beta-D-glucosiduronate is defined as the anion derived from beta-D-glucuronic acid.
    In beta-D-glucuronic acid the carboxylate group is located at the exocyclic CH2 (C6) 
    that is attached to a pyranose ring (sugar) and the aglycone attachment is via the 
    oxygen on the anomeric carbon (C1). 
    This function uses a SMARTS pattern that reflects the following key features:
      - A glycosidic –OC– linkage.
      - A pyranose ring with specific stereochemistry.
      - An exocyclic carboxylate (C(=O)[O-]) group.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a beta-D-glucosiduronate moiety, False otherwise.
        str: A reason for the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the beta-D-glucosiduronate moiety.
    # This pattern reflects:
    #   - "OC" at the beginning which is expected for the glycosidic linkage (the aglycone attached via an oxygen and a CH2 group)
    #   - A pyranose ring with defined beta stereochemistry.
    #   - An exocyclic carboxylate group attached to the ring (the characteristic feature of glucuronic acid).
    # Note: This SMARTS is heuristic; further refinement might be needed for edge cases.
    glucuronide_smarts = "OC[C@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@H]1C(=O)[O-]"
    glucuronide_pattern = Chem.MolFromSmarts(glucuronide_smarts)
    if glucuronide_pattern is None:
        return False, "Error in SMARTS pattern definition"
    
    # Check if the molecule has a substructure match for the beta-D-glucosiduronate moiety.
    if not mol.HasSubstructMatch(glucuronide_pattern):
        return False, "Beta-D-glucuronide moiety not found in the molecule"
    
    # Passed the substructure search.
    return True, "Contains beta-D-glucuronide moiety (anionic form of beta-D-glucuronic acid)"

# Example testing code (commented out):
# test_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@@]4([C@@H](CCC(O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C([O-])=O)=O)C)[H])[H])C)[H])C"
# result, reason = is_beta_D_glucosiduronate(test_smiles)
# print(result, reason)

# The provided test cases can be used to further validate the accuracy of this function.
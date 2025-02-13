"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group 
of any beta-D-glucuronic acid; major species at pH 7.3.
In glucuronic acid the primary alcohol group (CH2OH) of glucose is oxidized to a carboxylate 
group (C(=O)[O-]) and the aglycone is attached via the oxygen on the anomeric carbon.
This function uses a SMARTS pattern that (with relaxed stereochemistry) looks for:
  - A glycosidic oxygen linking the aglycone to the sugar.
  - A six-membered (pyranose) ring with at least three hydroxyl groups.
  - An exocyclic carboxylate group (C(=O)[O-]) attached to the ring.
Note that the SMARTS pattern is heuristic and may be refined further.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.

    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if a beta-D-glucosiduronate moiety is found, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the beta-D-glucosiduronate moiety.
    # The pattern (with no chiral specifications) is:
    #   O[C]1OC(O)C(O)C(O)C1C(=O)[O-]
    # which represents:
    #   - an oxygen (from the glycosidic bond)
    #   - a six-membered pyranose ring (atoms within 1...1) having one ring oxygen and four carbons substituted with OH groups
    #   - a carboxylate (C(=O)[O-]) attached to the ring (typically the exocyclic group from the original CH2OH).
    glucuronide_smarts = "O[C]1OC(O)C(O)C(O)C1C(=O)[O-]"
    glucuronide_pattern = Chem.MolFromSmarts(glucuronide_smarts)
    if glucuronide_pattern is None:
        return False, "Error in SMARTS pattern definition"
    
    # Check if the molecule has a matching substructure.
    if not mol.HasSubstructMatch(glucuronide_pattern):
        return False, "Beta-D-glucuronide moiety not found in the molecule"
    
    return True, "Contains beta-D-glucuronide moiety (anion of beta-D-glucuronic acid)"

# Example testing (commented out):
# test_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@@]4([C@@H](CCC(O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C([O-])=O)=O)C)[H])[H])C)[H])C"
# result, reason = is_beta_D_glucosiduronate(test_smiles)
# print(result, reason)
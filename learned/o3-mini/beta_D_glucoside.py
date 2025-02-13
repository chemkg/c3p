"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: Beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
This script uses a SMARTS substructure match to detect a beta-D-glucopyranoside fragment.
Note: Stereochemistry in sugars is subtle. The following SMARTS pattern is one approximation
that matches a common beta-D-glucoside ring fragment.
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    In this implementation the function looks for a beta-D-glucopyranoside substructure.
    The SMARTS pattern used is:
      $(OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O)
    which encodes a six-membered sugar ring (glucopyranose) with the glycosidic oxygen
    attached to its anomeric carbon in the beta configuration.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is identified as a beta-D-glucoside, False otherwise.
        str: Reason for the classification result.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for beta-D-glucopyranoside.
    # This pattern looks for a glucoside unit where the oxygen (glycosidic linker)
    # is bound to an anomeric carbon that has the proper stereochemistry (beta)
    beta_glucose_smarts = "[$(OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O)]"
    pattern = Chem.MolFromSmarts(beta_glucose_smarts)
    if pattern is None:
        return None, None

    # Search for the SMARTS pattern as a substructure in the molecule.
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return True, "Beta-D-glucoside substructure detected."
    else:
        return False, "No beta-D-glucoside fragment found in the molecule."
        
# Example usage (for testing purposes):
if __name__ == "__main__":
    # Testing with one of the examples: beta,beta-trehalose has two glycosidic bonds. 
    test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
    result, reason = is_beta_D_glucoside(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
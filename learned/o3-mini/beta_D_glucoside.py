"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: Beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
Improvements over the previous attempt:
  • Instead of using one very strict SMARTS string, we use several alternative
    patterns. Many beta-D-glucosides use either a direct O–glycosidic linkage or a 
    linkage via a CH2 (–O–C–) and the sugar may be in its six-membered (pyranose)
    or five-membered (furanose) form.
  • In the pyranose patterns we require the exocyclic CH2OH (characteristic of D-glucopyranose).
  • The chiral tags ([C@H] and [C@@H]) try to enforce the beta configuration at the anomeric center.
Note: Stereochemistry in sugars is very subtle so these SMARTS are only approximations.
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    This implementation searches for one of several alternative SMARTS patterns that
    each encode a beta-D-glucoside fragment. In particular:
      - For the pyranose (6-membered) form we look for a ring having an exocyclic CH2OH
        group and a beta-glycosidic oxygen attached either directly or via a CH2 linker.
      - For the furanose (5-membered) form we do the same.
    If any one of these substructures is found, the molecule is classified as a beta-D-glucoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a beta-D-glucoside fragment is detected, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a list of SMARTS patterns covering possible beta-D-glucoside fragments.
    # (a) Pyranose forms (6-membered sugar ring with an exocyclic CH2OH group).
    #     Pattern 1: direct attachment of the glycosidic oxygen.
    pattern_pyranose_direct = "[$(O[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1CO)]"
    #     Pattern 2: a CH2 linker between the aglycon oxygen and the sugar ring.
    pattern_pyranose_linker = "[$(OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1CO)]"
    
    # (b) Furanose forms (5-membered sugar ring with exocyclic CH2OH).
    #     Pattern 3: direct attachment.
    pattern_furanose_direct = "[$(O[C@H]1O[C@@H](O)[C@H](O)[C@H]1CO)]"
    #     Pattern 4: with a CH2 linker.
    pattern_furanose_linker = "[$(OC[C@H]1O[C@@H](O)[C@H](O)[C@H]1CO)]"
    
    # Compile the SMARTS patterns.
    smarts_patterns = []
    for pat_str in [pattern_pyranose_direct, pattern_pyranose_linker,
                    pattern_furanose_direct, pattern_furanose_linker]:
        pat = Chem.MolFromSmarts(pat_str)
        if pat is not None:
            smarts_patterns.append(pat)
    
    # If no pattern could be compiled, we cannot classify.
    if not smarts_patterns:
        return None, None
    
    # Check if any pattern matches the given molecule.
    for pat in smarts_patterns:
        if mol.HasSubstructMatch(pat):
            return True, "Beta-D-glucoside substructure detected."
    
    # If none of the patterns match, report a miss.
    return False, "No beta-D-glucoside fragment found in the molecule."
    
# Example usage for testing:
if __name__ == "__main__":
    # One test string from the provided examples.
    test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
    result, reason = is_beta_D_glucoside(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
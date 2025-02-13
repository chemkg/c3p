"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: pterocarpans
Definition:
  Members of the class of benzofurochromene with a 6a,11a-dihydro skeleton 
  (i.e. the 3,4-dihydro derivatives of coumestans) and its substituted derivatives.
  
  Our approach tries two SMARTS patterns:
    • Pattern 1 (pattern1) looks for an aromatic ring fused to two non‐aromatic (sp³) carbons,
      then an oxygen and a second aromatic ring. (This is similar to our first attempt.)
    • Pattern 2 (pattern2) relaxes the beginning of the SMARTS to allow the molecule to start with
      a tetrahedral (C;X4) center.
  If either pattern is found as a substructure, we classify the molecule as a (likely) pterocarpan.
  
  Note: This heuristic may not capture every pterocarpan and may give false negatives for some structures.
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule belongs to the pterocarpan class based on its SMILES string.
    The approach looks for a fused, tricyclic dihydrobenzofurochromene (pterocarpan) core.
    
    To maximize success we test for two similar substructure patterns:
      pattern1: "c1ccc2[C;!a][C;!a]Oc3ccccc3O2c1"
      pattern2: "[C;X4]1[C;X4]Oc2ccccc2O1"
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is (likely) a pterocarpan, else False.
        str : Explanation for the result.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the pterocarpan core.
    # Pattern 1 attempts to match a dihydrobenzofurochromene starting from an aromatic fragment.
    pattern_smarts1 = "c1ccc2[C;!a][C;!a]Oc3ccccc3O2c1"
    core1 = Chem.MolFromSmarts(pattern_smarts1)
    if core1 is None:
        return None, None  # Problem parsing SMARTS
    
    # Pattern 2 is an alternative, starting with a tetrahedral (sp3) carbon.
    pattern_smarts2 = "[C;X4]1[C;X4]Oc2ccccc2O1"
    core2 = Chem.MolFromSmarts(pattern_smarts2)
    if core2 is None:
        return None, None  # Problem parsing SMARTS
    
    # Check if molecule contains at least one of the core scaffolds.
    if mol.HasSubstructMatch(core1):
        return True, "Contains the pterocarpan core structure (matched pattern1: dihydrobenzofurochromene with aromatic starter)"
    if mol.HasSubstructMatch(core2):
        return True, "Contains the pterocarpan core structure (matched pattern2: alternative tetrahedral starting atom)"
    
    return False, "Pterocarpan core structure not found"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles_list = [
        # Test a few pterocarpan examples (stereochemistry removed for matching ease)
        "[H]C12COc3cc(O)c(CC=C(C)C)cc3C1Oc1c(CC=C(C)C)c(OC)c(O)cc21",  # lespeflorin G2
        "O1C2C(COC3=C2C=CC=4OC(C=CC34)(C)C)C5=C1C=C(O)C=C5",           # (-)-Shinpterocarpin
        "[H]C12COc3cc(O)ccc3C1Oc1cc3OCOc3cc21",                         # maackiain
        "COc1ccc2C3COc4cc(O)ccc4C3Oc2c1"                                  # medicarpin (simplified)
    ]
    for smi in test_smiles_list:
        res, reason = is_pterocarpans(smi)
        print(smi)
        print("Result:", res, "|", reason)
        print("-" * 60)
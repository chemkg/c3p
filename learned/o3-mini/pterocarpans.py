"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: pterocarpans
Definition:
  Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton 
  (i.e. the dihydro derivatives of coumestans) and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule belongs to the pterocarpan class based on its SMILES string.
    Pterocarpans are defined as having a benzofurochromene core with a 6a,11a-dihydro skeleton.
    The search is based on identifying a central scaffold that can be represented by the following
    dihydrofuran-fused bicyclic system:
         O1CCc2ccccc2OC1c3ccccc3
    This pattern is a simplified heuristic for a fused dihydrofuran ring (with one saturated carbon)
    bridging two aromatic rings. Many pterocarpan derivatives will have substituents on these rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is (likely) a pterocarpan, else False.
        str: Explanation for the classification result.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the pterocarpan core.
    # This pattern represents a dihydrofuran ring (O1CC...) fused to two benzene rings.
    # Note: This is a heuristic – not every pterocarpan may exactly match this query if heavily substituted.
    pattern_smarts = "O1CCc2ccccc2OC1c3ccccc3"
    core = Chem.MolFromSmarts(pattern_smarts)
    if core is None:
        return None, None  # In case the SMARTS pattern cannot be parsed
    
    # Check if the molecule contains the pterocarpan core scaffold.
    if not mol.HasSubstructMatch(core):
        return False, "Pterocarpan core structure not found"
    
    # (Optional additional checks can be done, such as molecular weight, rotatable bonds etc.
    #  but here we rely only on the substructure match to the core scaffold.)
    
    return True, "Contains the pterocarpan core structure (dihydrobenzofurochromene skeleton)"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "[H][C@]12COc3cc(O)ccc3[C@@]1([H])Oc1cc(OC)ccc21"  # (-)-medicarpin example
    result, reason = is_pterocarpans(test_smiles)
    print(result, reason)
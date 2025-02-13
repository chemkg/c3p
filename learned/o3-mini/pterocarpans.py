"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: pterocarpans
Definition:
  Members of the class of benzofurochromene with a 6a,11a-dihydro skeleton 
  (i.e. the 3,4-dihydro derivatives of coumestans) and its substituted derivatives.
  
  New approach: Instead of our previous simplified SMARTS, we now search for a core 
  pattern representing a fused tricyclic system where an aromatic ring (ring 1) is 
  fused to a second ring (central, containing two saturated carbons) that is in turn 
  fused (via an oxygen linkage) to a third aromatic ring. The SMARTS is a heuristic intended to 
  capture the dihydrobenzofurochromene framework observed in pterocarpans.
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule belongs to the pterocarpan class based on its SMILES string.
    Pterocarpans have a fused, tricyclic dihydrobenzofurochromene core. In our approach,
    we use a SMARTS query that looks for an aromatic ring (ring 1) fused to a central ring
    that contains two saturated (non-aromatic) carbons and an oxygen, and then fused to a
    second aromatic ring. This pattern is given below.

    SMARTS used: "c1ccc2[C;!a][C;!a]Oc3ccccc3O2c1"
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is (likely) a pterocarpan, else False.
        str : Explanation for the classification result.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern for the pterocarpan core.
    # This pattern looks for a fused tricyclic system:
    #  - "c1ccc2"           : an aromatic ring (ring 1) fused to ring 2.
    #  - "[C;!a][C;!a]"      : two saturated (non-aromatic) carbons (the dihydro centers, e.g. 6a and 11a).
    #  - "O"                 : an oxygen atom in the central ring.
    #  - "c3ccccc3"         : a fully aromatic ring (ring 3) fused to the central ring.
    #  - "O2" and "c1"       : closure of the fused rings.
    pattern_smarts = "c1ccc2[C;!a][C;!a]Oc3ccccc3O2c1"
    core = Chem.MolFromSmarts(pattern_smarts)
    if core is None:
        return None, None  # Unable to parse our SMARTS query
    
    # Check if the molecule contains the pterocarpan core scaffold.
    if not mol.HasSubstructMatch(core):
        return False, "Pterocarpan core structure not found"
    
    return True, "Contains the pterocarpan core structure (dihydrobenzofurochromene skeleton)"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: (-)-medicarpin simplified (remove stereochemistry for matching)
    test_smiles = "[H]C12COc3cc(O)ccc3C1Oc1cc(OC)ccc2"  
    result, reason = is_pterocarpans(test_smiles)
    print(result, reason)
"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: pterocarpans
Definition:
  Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton 
  (i.e. the 3,4-dihydro derivatives of coumestans) and its substituted derivatives.
  
  Our revised approach attempts to “lock in” the three fused rings. In this heuristic the pterocarpan core is 
  defined (after stripping off substituents) by two tetrahedral (sp³) carbons fused via an oxygen to an aromatic ring.
  
  We do this by requiring a match to a SMARTS pattern that is roughly equivalent to:
    [C;X4]1COC2c3ccccc3OC12
  This pattern requires that an sp³ carbon (atom with four connections) is attached to a CH₂ group which together 
  with an oxygen and an aromatic (benzene) ring complete two ring closures. (This is a coarse approximation but 
  it reduces the risk of false positives from an overly permissive alternative.)
  
  Note: The heuristic may miss out on certain distortions or unusual substitutions and is meant as a best‐effort filter.
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule belongs to the pterocarpan class based on its SMILES string.
    The approach uses a SMARTS pattern that attempts to capture a fused tricyclic core corresponding 
    approximately to the 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene scaffold.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is (likely) a pterocarpan, else False.
        str : Explanation for the result.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS that (heuristically) captures the pterocarpan fused tricyclic core.
    # This pattern requires:
    #   - a tetrahedral (sp³) carbon (via [C;X4]) at the start of a ring (ring closure label 1)
    #   - followed by a CH2 group (C) attached to an oxygen (O)
    #   - then a second ring (ring label 2) and an aromatic ring (c3ccccc3)
    #   - and finally the oxygens and ring closures to complete the tricyclic system.
    pattern_smarts = "[C;X4]1COC2c3ccccc3OC12"
    
    core = Chem.MolFromSmarts(pattern_smarts)
    if core is None:
        return None, "Error: Could not parse SMARTS pattern"
    
    # Check if the molecule matches the pterocarpan core pattern.
    if mol.HasSubstructMatch(core):
        return True, "Contains a fused tricyclic core consistent with a pterocarpan skeleton"
    
    # If no match, then it is likely not a pterocarpan
    return False, "Pterocarpan core structure not found"

# Example usage (for testing this module):
if __name__ == "__main__":
    test_smiles_list = [
        # Some examples from the provided list (stereochemistry simplified for matching)
        "[H][C@@]12COc3cc(O)c(CC=C(C)C)cc3[C@]1([H])Oc1c(CC=C(C)C)c(OC)c(O)cc21",  # lespeflorin G2
        "O1C2C(COC3=C2C=CC=4OC(C=CC34)(C)C)C5=C1C=C(O)C=C5",                       # (-)-Shinpterocarpin
        "[H]C12COc3cc(O)ccc3C1Oc1cc3OCOc3cc21",                                     # maackiain
        "O1C2C(C3=C1C=C(O)C=C3)COC4=C2C=C(C(O)=C4)CC=C(C)C",                         # Calocarpin
        # (Other examples from the list can be added here.)
    ]
    for smi in test_smiles_list:
        res, reason = is_pterocarpans(smi)
        print("SMILES:", smi)
        print("Result:", res, "|", reason)
        print("-" * 60)
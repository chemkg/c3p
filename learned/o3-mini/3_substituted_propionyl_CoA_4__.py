"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
Definition:
  An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
  of any 3-substituted propionyl-CoA; major species at pH 7.3.
  
Heuristic classification criteria:
  1. The molecule must be a valid structure.
  2. It must contain a thioester substructure (i.e. a carbonyl bound to sulfur, indicative of an acyl thioester).
  3. It must contain a characteristic CoA moiety fragment. Here we demand the presence of the pantetheine fragment
     “SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP” that is found in acyl-CoAs.
  4. It must have at least 4 deprotonated oxygen atoms ([O-]), matching its 4– overall charge.
  
If these conditions are satisfied, the molecule is taken to be a 3-substituted propionyl-CoA(4-).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule matches the 3-substituted propionyl-CoA(4-) criteria, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the presence of a thioester group.
    # SMARTS: a carbonyl ([C](=O)) directly attached to an S atom.
    thioester_smarts = "[C](=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester (C(=O)S) group found"
    
    # 2. Check for the CoA moiety.
    # Many acyl-CoAs contain the invariant pantetheine fragment:
    # "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP"
    # (This SMARTS is a simplified representation to capture a core part of the CoA scaffold.)
    coa_smarts = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety fragment not found"
    
    # 3. Check for the four deprotonated oxygens indicating the (4-) charge.
    # Count oxygen atoms (atomic number 8) with formal charge -1.
    deprot_oxy_count = sum(1 for atom in mol.GetAtoms() 
                             if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if deprot_oxy_count < 4:
        return False, f"Found {deprot_oxy_count} deprotonated oxygens; expected at least 4 for CoA(4-)"
    
    # If all checks pass, classify the molecule as a 3-substituted propionyl-CoA(4-)
    return True, "Molecule contains a thioester group, CoA moiety fragment, and sufficient deprotonated oxygens (4-)"
    
# Example usage:
if __name__ == '__main__':
    # You could test with one of the provided SMILES strings.
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    res, reason = is_3_substituted_propionyl_CoA_4__(test_smiles)
    print(res, reason)
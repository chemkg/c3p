"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
Definition:
  An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
  of any 3-substituted propionyl-CoA; major species at pH 7.3.

Heuristic classification criteria improved:
  1. The molecule must be a valid structure.
  2. It must contain a thioester group, i.e. a [C](=O)S fragment.
  3. It must contain a characteristic CoA moiety fragment (we search for the pantetheine core):
       SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP
  4. It must have at least 4 deprotonated oxygens ([O-]), consistent with a net (4-) charge.
  5. It must have exactly three carbonyl groups ([CX3](=O)): one thioester carbonyl plus the two carbonyls in the CoA fragment.
If all these conditions are met, we classify the molecule as a 3‑substituted propionyl-CoA(4-).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule conforms to the 3-substituted propionyl-CoA(4-) criteria.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the criteria, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a thioester group.
    # SMARTS: a carbonyl (C(=O)) directly attached to an S atom.
    thioester_smarts = "[C](=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester (C(=O)S) group found"
    
    # 2. Check for the CoA moiety by finding a core pantetheine fragment.
    # We use a simplified SMARTS capturing:
    # "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP"
    coa_smarts = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety fragment not found"
    
    # 3. Check for the four deprotonated oxygens ([O-]).
    deprot_oxy_count = sum(1 for atom in mol.GetAtoms() 
                             if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if deprot_oxy_count < 4:
        return False, f"Found {deprot_oxy_count} deprotonated oxygens; expected at least 4 for CoA(4-)"
    
    # 4. Count all carbonyl groups [CX3](=O) in the molecule.
    # In a proper 3-substituted propionyl-CoA(4-), there should be exactly three carbonyl groups:
    # one from the acyl thioester and two from the CoA fragment.
    carbonyl_smarts = "[CX3](=O)"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    n_carbonyl = len(carbonyl_matches)
    if n_carbonyl != 3:
        return False, f"Found {n_carbonyl} carbonyl groups; expected exactly 3 (1 acyl + 2 in CoA fragment)"
    
    return True, "Molecule contains a thioester group, CoA core fragment, sufficient deprotonated oxygens, and the expected number of carbonyl groups (3)"

# Example usage (for testing purposes):
if __name__ == '__main__':
    # Test with one of the provided SMILES strings (true positive example)
    test_smiles = "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    res, reason = is_3_substituted_propionyl_CoA_4__(test_smiles)
    print(res, reason)
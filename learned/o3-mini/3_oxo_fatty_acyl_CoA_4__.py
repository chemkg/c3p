"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI class “3-oxo-fatty acyl-CoA(4-)”
Definition: An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups 
of any 3-oxo-fatty acyl-CoA.

This version adds additional checks:
    • The overall formal charge must be –4.
    • The 3‑oxo fatty acyl thioester motif (SMARTS: "C(=O)CC(=O)S") is required.
      Also, the acyl (first) carbon in that motif must have exactly one carbonyl oxygen.
    • A CoA (pantetheine) signature is required (SMARTS: "SCCNC(=O)CCNC(=O)").
    • At least two phosphorus atoms exist.
    • Count of negatively charged oxygen atoms attached to phosphorus is at least 4.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise.
        str: A message providing the reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the overall formal charge is exactly –4.
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Overall formal charge is {total_charge}; expected –4 for a CoA(4–) derivative"
    
    # Look for the 3-oxo fatty acyl thioester motif.
    # The SMARTS "C(=O)CC(=O)S" is intended to capture an acyl carbonyl followed by a CH2 and then a second carbonyl attached to an S.
    oxo_smarts = "C(=O)CC(=O)S"
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "Does not contain the required 3-oxo fatty acyl thioester motif (C(=O)CC(=O)S)"
    
    # For each match, verify that the acyl (first) carbon has exactly one carbonyl oxygen.
    # (This check helps to rule out cases where extra carboxylate groups are present.)
    valid_oxo = False
    for match in oxo_matches:
        # In our SMARTS the first atom (index 0) is the acyl carbon.
        acyl_carbon = mol.GetAtomWithIdx(match[0])
        # Count oxygen neighbors connected by a double bond
        dO_count = 0
        for nbr in acyl_carbon.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(acyl_carbon.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                dO_count += 1
        if dO_count == 1:
            valid_oxo = True
            break
    if not valid_oxo:
        return False, "Acyl carbon in the 3-oxo motif does not have exactly one carbonyl oxygen (possible extra carboxylate group detected)"
    
    # Look for a typical CoA (pantetheine) moiety signature.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain the CoA moiety signature (SCCNC(=O)CCNC(=O))"
    
    # Check for the phosphorus atoms expected in a CoA derivative.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) < 2:
        return False, "Insufficient phosphorus atoms for a complete CoA moiety"
    
    # Count negatively charged oxygen atoms attached to phosphorus.
    neg_oxy_count = 0
    for p in phosphorus_atoms:
        for nbr in p.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == -1:
                neg_oxy_count += 1
    if neg_oxy_count < 4:
        return False, (f"Phosphate groups appear inadequately deprotonated (found {neg_oxy_count} negatively charged oxygens "
                       "attached to phosphorus; expected at least 4)")
    
    return True, ("Contains the 3-oxo fatty acyl thioester motif (with proper acyl carbon substitution), "
                  "a CoA moiety signature, and the expected deprotonated phosphate pattern (overall –4 charge)")

# For testing, one may uncomment and try a test SMILES:
# smiles_example = "CC(C)CCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
# result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
# print(result, reason)
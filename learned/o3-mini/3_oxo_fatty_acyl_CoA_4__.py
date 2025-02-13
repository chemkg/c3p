"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI class “3-oxo-fatty acyl-CoA(4-)”
Definition: An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups of any 3-oxo-fatty acyl-CoA.
This program uses SMARTS matching to look for:
    1. A 3-oxo fatty acyl chain motif (3‑oxo thioester): C(=O)CC(=O)S
    2. A CoA moiety signature from the pantetheine region: e.g., SCCNC(=O)CCNC(=O)
    3. Presence of deprotonated phosphate groups as indicated by oxygen atoms with a –1 formal charge 
       (we expect at least 4 such oxygens overall).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for the 3-oxo fatty acyl thioester motif using a SMARTS pattern.
    # This pattern looks for "C(=O)CC(=O)S".
    oxo_fatty_smarts = "C(=O)CC(=O)S"
    oxo_pattern = Chem.MolFromSmarts(oxo_fatty_smarts)
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Does not contain the required 3-oxo fatty acyl thioester motif (C(=O)CC(=O)S)"
    
    # Check for a typical CoA (pantetheine) signature.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain the CoA moiety signature (SCCNC(=O)CCNC(=O))"
    
    # Check for phosphorylation: rather than relying on the formal charge of phosphorus atoms,
    # we count oxygen atoms that show an explicit -1 formal charge and are bonded to phosphorus.
    # In a fully deprotonated CoA(4-) we expect at least 4 such oxygen atoms.
    neg_oxygen_count = 0
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) < 2:
        return False, "Insufficient phosphorus atoms for a complete CoA moiety"
    
    # For each phosphorus atom, count the oxygen neighbors with a formal charge of -1.
    for p in phosphorus_atoms:
        for nbr in p.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == -1:
                neg_oxygen_count += 1

    if neg_oxygen_count < 4:
        return False, ("Phosphate groups do not appear to be deprotonated adequately "
                       f"(found {neg_oxygen_count} negatively charged oxygen(s) attached to phosphorus, expected at least 4)")
    
    # Optional: check that the molecular weight is in an expected range for a CoA derivative.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for a CoA derivative"
    
    return True, ("Contains the 3-oxo fatty acyl chain motif, a characteristic CoA moiety signature, "
                  "and has the expected pattern of deprotonated phosphate groups")

# Example test (uncomment to perform a simple test)
# smiles_example = "CC(C)CCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
# result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
# print(result, reason)
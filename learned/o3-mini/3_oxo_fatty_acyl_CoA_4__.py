"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI class “3-oxo-fatty acyl-CoA(4-)”
Definition: An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups of any 
3-oxo-fatty acyl-CoA.
This program uses SMARTS matching to look for:
    1. A 3-oxo fatty acyl chain motif: a fragment matching C(=O)–C–C(=O)–S (i.e. a 3‑oxo thioester)
    2. A CoA moiety signature: a fragment typical of the pantetheine portion, e.g. SCCNC(=O)CCNC(=O)
    3. Presence of a sufficient number of phosphorus atoms that are deprotonated (formal charge negative)
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
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS for the 3-oxo fatty acyl motif.
    # This pattern seeks a substructure: carbonyl, then any carbon (usually CH2),
    # then another carbonyl attached to a sulfur.
    oxo_fatty_smarts = "C(=O)CC(=O)S"
    oxo_pattern = Chem.MolFromSmarts(oxo_fatty_smarts)
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Does not contain the required 3-oxo fatty acyl thioester motif (C(=O)CC(=O)S)"
    
    # Define a SMARTS pattern to catch a CoA signature from the pantetheine region.
    # This pattern looks for a fragment: S - C - C - N - C(=O) - C - C - N - C(=O)
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain the CoA moiety signature (SCCNC(=O)CCNC(=O))"
    
    # Count phosphorus atoms that are likely present in the phosphate groups.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    # Among phosphorus, check if at least one has a negative formal charge.
    neg_P_count = sum(1 for atom in phosphorus_atoms if atom.GetFormalCharge() < 0)
    if len(phosphorus_atoms) < 2:
        return False, "Insufficient phosphorus atoms for a complete CoA moiety"
    if neg_P_count < 2:
        return False, "Phosphate groups do not appear to be deprotonated (not enough negatively charged phosphorus atoms)"
    
    # Additional optional check: molecular weight might be in an expected range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight is lower than expected for a CoA derivative"
    
    return True, "Contains the 3-oxo-fatty acyl chain motif, a CoA moiety signature, and appropriate phosphate groups"

# Example test (uncomment to run a simple test)
# smiles_example = "CC(C)CCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
# result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
# print(result, reason)
"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol (a glycerophosphoglycerol)
Definition: “A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.”
This program checks for the presence of (1) phosphorus, (2) at least two ester groups 
(as a proxy for two acyl chains) and (3) a glycerol phosphate head group substructure.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol (PG)
    based on its SMILES string.
    
    A phosphatidylglycerol should have:
      - A valid chemical structure.
      - At least one phosphorus atom.
      - At least two acyl (ester) groups (OC(=O)) indicating diacyl chains.
      - A glycerol phosphate head-group (an –O–P(=O)(O)OC(CO)O fragment).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple containing True with a success message if the molecule 
                     meets the criteria; False with a reason otherwise.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of phosphorus element
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found; cannot be a phosphatidylglycerol"
    
    # Count ester groups (this simple pattern “OC(=O)” serves as a proxy for acyl chain attachments)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester group(s); need at least 2 for diacyl chains"
    
    # Check for a glycerol phosphate head-group.
    # This SMARTS looks for an oxygen bonded to a phosphorus that is double bonded to oxygen and connected through oxygen to a glycerol-like fragment.
    # Here, we ignore stereochemical details.
    headgroup_pattern = Chem.MolFromSmarts("O-P(=O)(O)OC(CO)O")
    if not mol.HasSubstructMatch(headgroup_pattern):
        return False, "Glycerol phosphate head group not found"
    
    # Optionally, one can check for a minimum molecular weight or long aliphatic chains, but many PGs are diverse.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # heuristic cutoff
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical phosphatidylglycerol"
    
    return True, "Molecule contains a glycerol phosphate head group with diacyl chains, consistent with phosphatidylglycerol"

# Example usage (uncomment for testing):
# smiles_example = "C([C@@](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O"
# result, reason = is_phosphatidylglycerol(smiles_example)
# print(result, reason)
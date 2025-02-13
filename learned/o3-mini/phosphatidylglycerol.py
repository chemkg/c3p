"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol (a glycerophosphoglycerol)
Definition: “A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.”
This program attempts to identify a bona fide phosphatidylglycerol by ensuring:
  1) The molecule is valid.
  2) It contains exactly one phosphorus (typical of a single PG head group).
  3) It contains exactly two ester groups (as a proxy for two acyl chains attached at sn-1 and sn-2).
  4) It contains a glycerol phosphate headgroup linked to a free glycerol, here captured by the SMARTS:
       "OP(=O)(O)OC(CO)O"
     (this pattern tries to match a phosphate group attached via oxygen to a glycerol moiety).
  5) Its molecular weight is above a heuristic cutoff.
Note: PG structures are diverse. No simple set of substructure rules is perfect, but these
criteria help to reduce mis‐classifications.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol (PG) based on its SMILES string.
    
    A phosphatidylglycerol should have:
      - A valid chemical structure.
      - Exactly one phosphorus atom.
      - Exactly two ester groups (OC(=O)) indicating diacyl chains.
      - A glycerol phosphate head-group, in which a phosphate is linked to a glycerol 
        that in turn is linked via an –O– bond to another glycerol. 
        (We use the SMARTS "OP(=O)(O)OC(CO)O" as a proxy.)
      - A heuristic minimum molecular weight.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        (bool, str): Tuple; True with a success message if criteria are met,
                     False with an explanation otherwise.
    """
    # Parse SMILES; if invalid, return early.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that exactly one phosphorus is present.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly 1 phosphorus atom; found {len(p_atoms)}"
    
    # Count ester groups (using a SMARTS for OC(=O) as a proxy for acyl chains).
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Expected exactly 2 ester groups; found {len(ester_matches)}"
    
    # Check for the glycerol phosphate head-group.
    # This pattern is meant to capture a phosphate (OP(=O)(O)) connected via oxygen
    # to a glycerol fragment (OC(CO)O). (Stereo details are ignored.)
    headgroup_pattern = Chem.MolFromSmarts("OP(=O)(O)OC(CO)O")
    headgroup_matches = mol.GetSubstructMatches(headgroup_pattern)
    if len(headgroup_matches) != 1:
        return False, "Glycerol phosphate head-group not found exactly once"
    
    # Optionally, check for a minimum molecular weight for phosphatidylglycerol structures.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # heuristic cutoff
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for typical phosphatidylglycerol"
    
    return True, "Molecule has a single phosphorus, 2 acyl (ester) chains, and a glycerol phosphate head-group consistent with phosphatidylglycerol"

# Example usage (for testing purposes):
# smiles_example = "C([C@@](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O"
# result, reason = is_phosphatidylglycerol(smiles_example)
# print(result, reason)
"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol (a glycerophosphoglycerol)
Definition: “A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.”
This classifier attempts to identify bona fide phosphatidylglycerols by enforcing:
  1) A valid molecule.
  2) Exactly one phosphorus atom.
  3) Exactly two acyl ester groups (matching the SMARTS "OC(=O)").
  4) A glycerol phosphate head‐group, here approximated by the SMARTS "P(=O)(O)OCC(CO)O"
     (which should capture a phosphate attached via an oxygen to a glycerol fragment).
  5) A molecular weight in a reasonable range for typical PG molecules (600–1000 Da).
Note: Many PG subclasses exist, and this set of rules is heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol (PG) based on its SMILES string.
    
    PGs are defined here as diacylglycerols whose free hydroxyl substituent has been replaced
    by a phosphatidyl group (i.e. a phosphate attached to a glycerol head-group). Thus we require:
      - A valid chemical structure.
      - Exactly one phosphorus atom.
      - Exactly two acyl ester groups (as a proxy for the two fatty acid chains).
      - A glycerol phosphate head‐group with the connectivity P(=O)(O)OCC(CO)O (ignoring stereochemistry).
      - A molecular weight within a heuristic range for PG (600 –1000 Da).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple; True with a success message if criteria are met,
                     otherwise False with an explanation.
    """
    # Parse SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one phosphorus atom.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly 1 phosphorus atom; found {len(p_atoms)}"
    
    # Count acyl ester groups using the SMARTS "OC(=O)".
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Expected exactly 2 acyl ester groups; found {len(ester_matches)}"
    
    # Check for glycerol phosphate head-group.
    # We use a SMARTS pattern designed to match a phosphate bound via an oxygen to
    # a glycerol-like moiety. We ignore chirality.
    headgroup_pattern = Chem.MolFromSmarts("P(=O)(O)OCC(CO)O")
    headgroup_matches = mol.GetSubstructMatches(headgroup_pattern, useChirality=False)
    if len(headgroup_matches) < 1:
        return False, "Glycerol phosphate head-group not found"
    if len(headgroup_matches) > 1:
        return False, "Multiple glycerol phosphate head-group patterns found"
    
    # Check molecular weight is in a typical range for PG molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (600 <= mol_wt <= 1000):
        return False, f"Molecular weight {mol_wt:.1f} Da out of typical PG range (600–1000 Da)"
    
    return True, ("Molecule has a single phosphorus, 2 acyl (ester) chains, and a glycerol phosphate head-group "
                  "consistent with phosphatidylglycerol")

# Example usage for testing (uncomment to run):
# smiles_example = "C([C@@](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O"
# result, reason = is_phosphatidylglycerol(smiles_example)
# print(result, reason)
"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol (a glycerophosphoglycerol)
Definition: “A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.”
Heuristic rules used:
  1) Must be a valid molecule.
  2) Should have one or two phosphorus atoms (PGP typically shows two P, while PG shows one).
     Molecules with >2 phosphorus atoms (e.g. cardiolipin) are rejected.
  3) Should have exactly two acyl ester groups (using the SMARTS "OC(=O)" to represent the ester bonds).
  4) Must contain a glycerol phosphate head‐group.
     We use a SMARTS pattern (ignoring chirality) designed to capture a phosphate 
     bound via an oxygen to a glycerol fragment.
  5) No nitrogen atoms are allowed.
  6) Molecular weight typical for PG (600–1000 Da).
Note: This is a heuristic and may miss or wrongly classify molecules in sidebar cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol (PG) or a closely related PGP based on its SMILES string.
    
    PGs are defined here as diacylglycerols whose free hydroxyl substituent has been replaced
    by a phosphatidyl group (i.e. a phosphate attached to a glycerol head-group). Thus we require:
      - A valid chemical structure.
      - 1 or 2 phosphorus atoms (reject if phosphorus count is 0 or >2).
      - Exactly two acyl ester groups (using the SMARTS "OC(=O)" as a proxy for the fatty acyl attachments).
      - A glycerol phosphate head‐group; we look for a fragment resembling "OP(=O)(O)OCC(O)CO" (ignoring chirality).
      - No nitrogen atoms (PGs do not contain heteroatoms like N).
      - A molecular weight in a typical PG range (600 – 1000 Da).
      
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
    
    # Check molecule contains no nitrogen atoms.
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if n_atoms:
        return False, "Contains nitrogen atoms, not typical for PG"
    
    # Check for phosphorus atoms: allow 1 or 2 (PG vs. PGP)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not (1 <= len(p_atoms) <= 2):
        return False, f"Expected 1 or 2 phosphorus atoms; found {len(p_atoms)}"
    
    # Count acyl ester groups. These bonds are usually represented by the pattern "OC(=O)"
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Expected exactly 2 acyl ester groups; found {len(ester_matches)}"
    
    # Look for a glycerol phosphate head-group.
    # The pattern below is designed (ignoring chirality) to capture a phosphate bound via an oxygen
    # to a glycerol-like fragment. We do not insist on exactly one match since overlapping substructures
    # can sometimes give multiple hits.
    headgroup_pattern = Chem.MolFromSmarts("OP(=O)(O)OCC(O)CO")
    headgroup_matches = mol.GetSubstructMatches(headgroup_pattern, useChirality=False)
    if len(headgroup_matches) < 1:
        return False, "Glycerol phosphate head-group not found"
    
    # Check the molecular weight is within typical PG range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (600 <= mol_wt <= 1000):
        return False, f"Molecular weight {mol_wt:.1f} Da out of typical PG range (600–1000 Da)"
    
    return True, ("Molecule has the expected phosphorus count (1–2), 2 acyl ester chains, a glycerol phosphate head-group, "
                  "and no nitrogen atoms with molecular weight in the PG range")
    
# Example usage for testing (uncomment to run):
# smiles_example = "C([C@@](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O"
# result, reason = is_phosphatidylglycerol(smiles_example)
# print(result, reason)
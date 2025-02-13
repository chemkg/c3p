"""
Classifies: CHEBI:17334 penicillin
"""
# penicillin.py
"""
Classifies: Penicillin class
Definition:
  Any member of the group of substituted penams containing two methyl substituents at position 2,
  a carboxylate substituent at position 3 and a carboxamido group at position 6.
  
This heuristic looks for:
  • A beta‐lactam substructure (four‐membered cyclic amide)
  • A [S]-C fragment with two methyl substituents on the C (indicative of C2 dimethyl substitution)
  • A carboxylate group (–C(=O)O or –C(=O)[O–])
  • A carboxamido group (–NC(=O)...)
  
Note: The penicillin scaffold is quite complex; the following code uses several SMARTS patterns
as a heuristic for a “penicillin” structure.
"""

from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    
    The function uses SMARTS patterns to look for:
      1. A beta-lactam ring (a 4-membered ring with a cyclic amide)
      2. A sulfur-bound carbon having two methyl groups (indicative of the 2,2-dimethyl substituents)
      3. A carboxylate substituent (C(=O)O or C(=O)[O-])
      4. A carboxamido group (NC(=O)...)
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as penicillin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a beta-lactam ring.
    # A simple SMARTS for a 4-membered cyclic amide: ring of N-C-C-C with a carbonyl.
    beta_lactam_smarts = "N1C(=O)CC1"
    beta_lactam_query = Chem.MolFromSmarts(beta_lactam_smarts)
    if not mol.HasSubstructMatch(beta_lactam_query):
        return False, "Missing beta-lactam (4-membered cyclic amide) ring"
    
    # 2. Check for a two-methyl pattern at position 2.
    # We look for an S-atom bonded to a carbon that carries two methyl groups.
    # This SMARTS looks for an S attached to a carbon that in turn is bound to two CH3 groups.
    methyl_at_pos2_smarts = "[S]-C([CH3])([CH3])"
    methyl_query = Chem.MolFromSmarts(methyl_at_pos2_smarts)
    if not mol.HasSubstructMatch(methyl_query):
        return False, "Missing S-bound carbon with two methyl substituents (expected at position 2)"
    
    # 3. Check for a carboxylate substituent.
    # We allow either protonated carboxylic acid or deprotonated form. 
    carboxylate_smarts1 = "C(=O)[O-]"
    carboxylate_smarts2 = "C(=O)O"
    carboxylate_query1 = Chem.MolFromSmarts(carboxylate_smarts1)
    carboxylate_query2 = Chem.MolFromSmarts(carboxylate_smarts2)
    if not (mol.HasSubstructMatch(carboxylate_query1) or mol.HasSubstructMatch(carboxylate_query2)):
        return False, "Missing carboxylate (C(=O)O or C(=O)[O-]) substituent (expected at position 3)"
    
    # 4. Check for a carboxamido group.
    # We look for an amide where the nitrogen is directly bound (i.e. NC(=O)).
    carboxamido_smarts = "NC(=O)"
    carboxamido_query = Chem.MolFromSmarts(carboxamido_smarts)
    if not mol.HasSubstructMatch(carboxamido_query):
        return False, "Missing carboxamido group (NC(=O)-, expected at position 6)"
    
    return True, "Matches penicillin criteria: Contains beta‐lactam, 2,2-dimethyl, carboxylate and carboxamido substituents"

# For testing purposes, one could call the function with one of the example SMILES strings.
# For example:
# result, reason = is_penicillin("[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)CCC[C@H](N)C(O)=O)")
# print(result, reason)
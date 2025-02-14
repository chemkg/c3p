"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for CoA core structure, including the deprotonated phosphates
    #  Adenosine diphosphate part (minus one phosphate): [C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP([O-])(=O)OP([O-])=O)[C@H]([C@H]1O)
    #  Pantetheine part: NCC(=O)CCNC(=O)[C@H](O)C(C)(C)
    #  mercaptoethylamine: -S-C-C-N 
    # relaxed core substructure, with ring carbon substitution on ribose
    coa_core_pattern = Chem.MolFromSmarts("[NX2]c1nc(N)nc[n]1[C@H]1[C@@H]([C@@H]([C@H]([CH]1)OP([O-])(=O)[OP]([O-])([O-])=O)O)O")
    if not mol.HasSubstructMatch(coa_core_pattern):
        return False, "CoA core structure not found"

    #  Check for thioester bond (-S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("[#16][CX3](=[OX1])")
    if not mol.HasSubstructMatch(thioester_pattern):
      return False, "Thioester bond not found"


    # 2. Check for long-chain fatty acid part (more than 10 carbons, could be branched or not)
    # This is more complex because the chain can be variable. We need at least 10 carbons.
    # We also consider double bonds as part of the chain.
    
    # Look for a long carbon chain with at least 10 carbons with possible double bonds
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
         return False, "No long carbon chain found in the fatty acid part"
    
    # Count the carbons in the fatty acyl chain
    fatty_acid_carbon_count = 0
    # Find the thioester bond
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester found - cannot find fatty acid."
    
    thioester_start = thioester_matches[0][1]
    for bond in mol.GetAtomWithIdx(thioester_start).GetBonds():
        if bond.GetBeginAtomIdx() == thioester_start:
            next_atom = bond.GetEndAtomIdx()
        else:
            next_atom = bond.GetBeginAtomIdx()
        if mol.GetAtomWithIdx(next_atom).GetAtomicNum() == 6:
             fatty_acid_start_atom = next_atom
             break
    else:
        return False, "Could not trace fatty acid chain"
        
    # Traverse the chain to count carbons
    visited = {fatty_acid_start_atom}
    q = [fatty_acid_start_atom]
    
    while q:
        current_atom_idx = q.pop(0)
        
        if mol.GetAtomWithIdx(current_atom_idx).GetAtomicNum() == 6:
              fatty_acid_carbon_count +=1
            
        for bond in mol.GetAtomWithIdx(current_atom_idx).GetBonds():
          if bond.GetBeginAtomIdx() == current_atom_idx:
            next_atom = bond.GetEndAtomIdx()
          else:
             next_atom = bond.GetBeginAtomIdx()
          
          if next_atom not in visited and mol.GetAtomWithIdx(next_atom).GetAtomicNum() in [6,1]: #Count also hydrogens in chain
             q.append(next_atom)
             visited.add(next_atom)
    
    if fatty_acid_carbon_count < 10 :
        return False, f"Fatty acid chain is too short {fatty_acid_carbon_count} carbons."


    # 3. Check for 4- charge (three deprotonated phosphates)
    # Count number of deprotonated O atoms (charge -1) in phosphates
   # Note, for these molecules we need to check if there are 3 phosphate groups deprotonated
    phosphate_pattern = Chem.MolFromSmarts("[P]([O-])(=O)([O-])")
    num_phosphates = len(mol.GetSubstructMatches(phosphate_pattern))
    if num_phosphates < 2:
        return False, "Not enough deprotonated phosphates"
    
    mol_charge = rdMolDescriptors.CalcMolCharge(mol)
    if mol_charge != -4:
        return False, f"Molecule does not have -4 charge, found {mol_charge}"

    return True, "Molecule is a long-chain fatty acyl-CoA(4-)."
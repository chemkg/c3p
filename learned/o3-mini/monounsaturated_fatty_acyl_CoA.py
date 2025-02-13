"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated Fatty Acyl-CoA
Definition: An unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon‐carbon double bond.
Our strategy is:
  1. Verify the molecule has a CoA moiety (by requiring the adenine group) 
     and a thioester linkage (C(=O)S) connecting a fatty acyl chain to the CoA.
  2. Identify the thioester bond and then “disconnect” the fatty acyl chain by breaking
     the bond that links the carbonyl to its R‐group.
  3. Extract that fragment (the fatty acyl chain) and count its carbon–carbon double bonds.
     The acyl chain must have exactly one C=C bond.
"""

from rdkit import Chem
from rdkit.Chem import rdchem, rdmolops

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    The method:
      1. Checks for a CoA signature (via an adenine substructure) and a thioester bond.
      2. Finds the thioester bond (C(=O)S) and determines (via the carbonyl carbon)
         which substituent is the acyl chain.
      3. “Disconnects” the fatty acyl chain from the CoA by breaking that bond,
         extracts the fragment corresponding to the fatty acyl chain,
         and counts the number of C=C bonds (between carbon atoms).
      4. Returns True only if exactly one C=C bond is found in that fragment.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise.
        str: Reason for classification or error.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for CoA signature.
    # We require an adenine ring, which is present in almost every CoA derivative.
    adenine_smarts = "N1C=NC2=C1N=CN=C2"  # a common adenine fragment
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine_query):
        return False, "CoA moiety (adenine substructure) not detected"
    
    # Step 2: Look for a thioester bond using SMARTS "C(=O)S"
    thioester_smarts = "C(=O)S"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_query)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond found; not a fatty acyl-CoA"
    
    # Assume the first thioester match corresponds to the fatty acyl chain linkage.
    # In the SMARTS "C(=O)S", match[0] is the carbonyl carbon and match[1] is the sulfur.
    carbonyl_idx = thioester_matches[0][0]
    sulfur_idx    = thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the neighbor of the carbonyl atom that is part of the fatty acyl chain.
    # The carbonyl will have (typically) two neighbors: one oxygen from C=O and one group (the acyl chain).
    acyl_start_idx = None
    for nbr in carbonyl_atom.GetNeighbors():
        # Exclude the oxygen (its atomic number is 8)
        if nbr.GetAtomicNum() == 8:
            continue
        # Exclude the sulfur atom belonging to the thioester
        if nbr.GetIdx() == sulfur_idx:
            continue
        # The remaining neighbor is assumed to be the start of the fatty acyl chain.
        acyl_start_idx = nbr.GetIdx()
        break
    if acyl_start_idx is None:
        return False, "No fatty acyl chain attached to the thioester carbonyl found"
    
    # Now we want to 'disconnect' the fatty acyl chain.
    # We do this by breaking the bond between the carbonyl and the acyl_start atom.
    # Find the bond between carbonyl atom and acyl_start:
    bond_to_break = mol.GetBondBetweenAtoms(carbonyl_idx, acyl_start_idx)
    if bond_to_break is None:
        return False, "Bond between carbonyl and fatty acyl chain not found"
    bond_idx = bond_to_break.GetIdx()
    
    # Break this bond in the molecule, introducing dummy atoms.
    # The new molecule (frag_mol) will have two additional dummy atoms at the break sites.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    
    # Now, get the fragments as tuples of the original atom indices.
    # The fragment corresponding to the fatty acyl chain is the one that contains acyl_start_idx.
    frags = rdmolops.GetMolFrags(frag_mol, asMols=False, sanitizeFrags=False)
    acyl_frag_indices = None
    for frag in frags:
        if acyl_start_idx in frag:
            acyl_frag_indices = frag
            break
    if acyl_frag_indices is None:
        return False, "Could not isolate the fatty acyl chain fragment"
    
    # Create a new molecule from the acyl fragment.
    acyl_chain_mol = Chem.PathToSubmol(mol, acyl_frag_indices)
    
    # For debugging one might want to see the SMILES for the acyl fragment:
    # acyl_smiles = Chem.MolToSmiles(acyl_chain_mol)
    # print("Acyl chain fragment:", acyl_smiles)
    
    # Count the number of carbon-carbon double bonds in the acyl chain fragment.
    double_bond_count = 0
    for bond in acyl_chain_mol.GetBonds():
        # Consider only bonds connecting two carbon atoms.
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                double_bond_count += 1
    
    # According to our definition the fatty acyl chain must have exactly one C=C bond.
    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} C=C bond(s) instead of exactly one"
    
    return True, "Molecule is a monounsaturated fatty acyl-CoA with one C=C bond in the fatty acyl chain"


# For testing purposes:
if __name__ == '__main__':
    test_smiles = "CCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    valid, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(valid, "->", reason)
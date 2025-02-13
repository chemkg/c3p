"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated Fatty Acyl-CoA
Definition: Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.
Strategy:
  1. Check that the molecule contains a CoA moiety by matching an adenine substructure (using an aromatic SMARTS pattern)
     and that it has a thioester linkage (C(=O)S) that connects the fatty acyl chain to the CoA.
  2. Locate the thioester bond and from its carbonyl atom isolate the fatty acyl chain.
  3. "Disconnect" the acyl chain by breaking the bond between the carbonyl and its alkyl substituent.
  4. Count the number of C=C bonds in the acyl chain fragment. The chain must have exactly one such bond.
"""

from rdkit import Chem
from rdkit.Chem import rdchem, rdmolops

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    The classification uses the following steps:
      1. Verifies that the molecule has a CoA signature by checking for an adenine substructure
         (using aromatic SMARTS) and a thioester bond (C(=O)S).
      2. Finds the thioester bond and identifies the fatty acyl chain attached to the carbonyl.
      3. "Disconnects" the acyl chain from the rest of the molecule.
      4. Counts the carbon-carbon double bonds in the isolated acyl chain; the molecule qualifies only if exactly one is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for the CoA signature.
    # Use an aromatic SMARTS pattern for adenine (lower-case letters ensure matching aromatic nitrogens).
    adenine_smarts = "n1cnc2c(n1)ncnc2"
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine_query):
        return False, "CoA moiety (adenine substructure) not detected"
    
    # Check for the thioester bond (SMARTS for C(=O)S).
    thioester_smarts = "C(=O)S"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_query)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond found; not a fatty acyl-CoA"
    
    # Assume that the first thioester match corresponds to the fatty acyl chain linkage.
    # In the match, index 0 is the carbonyl carbon, index 1 is the attached sulfur.
    carbonyl_idx = thioester_matches[0][0]
    sulfur_idx    = thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the neighbor of the carbonyl atom that is part of the fatty acyl chain.
    # Exclude the oxygen from the carbonyl (atomic number 8) and the sulfur from the thioester.
    acyl_start_idx = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 8:  # skip oxygen
            continue
        if nbr.GetIdx() == sulfur_idx:  # skip sulfur
            continue
        # Found the neighbor that acts as the start of the fatty acyl chain.
        acyl_start_idx = nbr.GetIdx()
        break
    if acyl_start_idx is None:
        return False, "No fatty acyl chain attached to the thioester carbonyl found"
    
    # Step 2: Disconnect the fatty acyl chain.
    # Break the bond between the carbonyl carbon and the acyl chain start atom.
    bond_to_break = mol.GetBondBetweenAtoms(carbonyl_idx, acyl_start_idx)
    if bond_to_break is None:
        return False, "Bond between carbonyl and fatty acyl chain not found"
    bond_idx = bond_to_break.GetIdx()
    
    # Fragment the molecule at the identified bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    
    # Get fragments (each as a tuple of original atom indices) and locate the fragment containing the acyl chain start.
    frags = rdmolops.GetMolFrags(frag_mol, asMols=False, sanitizeFrags=False)
    acyl_frag_indices = None
    for frag in frags:
        if acyl_start_idx in frag:
            acyl_frag_indices = frag
            break
    if acyl_frag_indices is None:
        return False, "Could not isolate the fatty acyl chain fragment"
    
    # Create a new molecule from the isolated fatty acyl chain fragment.
    acyl_chain_mol = Chem.PathToSubmol(mol, acyl_frag_indices)
    
    # Step 3: Count the number of carbon-carbon double bonds in the acyl chain fragment.
    double_bond_count = 0
    for bond in acyl_chain_mol.GetBonds():
        # Only consider bonds between carbon atoms.
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                double_bond_count += 1
                
    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} C=C bond(s) instead of exactly one"
    
    return True, "Molecule is a monounsaturated fatty acyl-CoA with one C=C bond in the fatty acyl chain"

# For testing purposes:
if __name__ == '__main__':
    test_smiles = "CCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    valid, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(valid, "->", reason)
"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated Fatty Acyl-CoA
Definition: Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.
Strategy:
  1. Verify that the molecule contains a CoA moiety by checking for an adenine substructure.
     Instead of relying on one SMARTS pattern, use a couple of alternatives.
  2. Check that a thioester bond (C(=O)S) is present.
  3. Identify the thioester bond, then find and "disconnect" the fatty acyl chain attached via the carbonyl.
  4. Count the number of carbon-carbon double bonds in the resulting acyl chain fragment.
     Accept the molecule only if exactly one C=C bond is found.
"""

from rdkit import Chem
from rdkit.Chem import rdchem, rdmolops

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    
    The classification uses the following steps:
      1. Parses the SMILES string.
      2. Checks for the presence of a CoA signature by looking for an adenine substructure.
         (Several SMARTS for adenine are tried.)
      3. Checks for the thioester bond (C(=O)S) linking the fatty acyl portion.
      4. Disconnects the fatty acyl chain from the rest of the molecule.
      5. Counts the carbon-carbon double bonds in the isolated acyl chain; exactly one is required.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for the CoA signature (adenine substructure).
    # Try several SMARTS representations to capture variations in adenine.
    adenine_smarts_patterns = ["n1cnc2c(n1)ncnc2", "N1C=NC2=C1N=CN=C2N", "c1nc2c(cn1)nc(n2)"]
    adenine_found = False
    for smarts in adenine_smarts_patterns:
        query = Chem.MolFromSmarts(smarts)
        if query and mol.HasSubstructMatch(query):
            adenine_found = True
            break
    if not adenine_found:
        return False, "CoA moiety (adenine substructure) not detected"
    
    # Step 2: Check for the thioester bond: SMARTS for C(=O)S.
    thioester_smarts = "C(=O)S"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_query)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond found; not a fatty acyl-CoA"
    
    # Assume that the first thioester match corresponds to the fatty acyl linkage.
    # In this match, index 0 is the carbonyl carbon, index 1 is the attached sulfur.
    carbonyl_idx = thioester_matches[0][0]
    sulfur_idx    = thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Step 3: Identify the neighbor of the carbonyl atom that is part of the fatty acyl chain.
    # Skip the oxygen of the carbonyl and the sulfur of the thioester.
    acyl_start_idx = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 8:  # skip oxygen
            continue
        if nbr.GetIdx() == sulfur_idx:  # skip sulfur
            continue
        acyl_start_idx = nbr.GetIdx()
        break
    if acyl_start_idx is None:
        return False, "No fatty acyl chain attached to thioester carbonyl found"
    
    # Step 4: Disconnect the fatty acyl chain.
    bond_to_break = mol.GetBondBetweenAtoms(carbonyl_idx, acyl_start_idx)
    if bond_to_break is None:
        return False, "Bond between carbonyl and fatty acyl chain not found"
    bond_idx = bond_to_break.GetIdx()
    
    # Fragment the molecule at the identified bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    
    # Get fragments (each as a tuple of atom indices) and locate the fragment containing the acyl chain start.
    frags = rdmolops.GetMolFrags(frag_mol, asMols=False, sanitizeFrags=False)
    acyl_frag_indices = None
    for frag in frags:
        if acyl_start_idx in frag:
            acyl_frag_indices = frag
            break
    if acyl_frag_indices is None:
        return False, "Could not isolate the fatty acyl chain fragment"
    
    # Use PathToSubmol to extract the fragment as a separate molecule.
    acyl_chain_mol = Chem.PathToSubmol(mol, acyl_frag_indices)
    
    # Step 5: Count the number of carbon-carbon double bonds in the acyl chain fragment.
    double_bond_count = 0
    for bond in acyl_chain_mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # consider only bonds between carbon atoms
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                double_bond_count += 1
                
    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} C=C bond(s) instead of exactly one"
    
    return True, "Molecule is a monounsaturated fatty acyl-CoA with one C=C bond in the fatty acyl chain"

# For testing purposes:
if __name__ == '__main__':
    # Example test: (11Z)-hexadec-11-enoyl-CoA
    test_smiles = "CCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    valid, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(valid, "->", reason)
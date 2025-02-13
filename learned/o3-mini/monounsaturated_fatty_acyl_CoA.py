"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated Fatty Acyl-CoA
Definition: An unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.
This function looks for a thioester bond (C(=O)S) linking a fatty acyl chain to the remainder of the molecule (the CoA part).
Then it uses a breadth-first search (excluding the carbonyl carbon) to collect the atoms belonging to the fatty acyl chain and counts
the number of carbon-carbon double bonds within that fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    In our strategy:
      1. The molecule must have a thioester bond (i.e. a C(=O)S substructure) that connects the fatty acyl chain to CoA.
      2. One of the substituents of the carbonyl (C in C(=O)) must be the fatty acyl chain. We then traverse from that atom,
         deliberately not backtracking into the carbonyl so as to “isolate” the acyl chain.
      3. The acyl chain (the carbons reached in that search) should have exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule has a thioester bond (C(=O)S).
    thioester_smarts = "C(=O)S"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_query)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond found; not a fatty acyl-CoA"

    # For simplicity we assume the first match is the fatty acyl-CoA linkage.
    # In the pattern "C(=O)S", match[0] is the carbonyl carbon and match[1] is the thioester sulfur.
    carbonyl_idx = thioester_matches[0][0]
    sulfur_idx   = thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)

    # In a fatty acyl group, the carbonyl carbon has two substituents: one oxygen from the C=O and one group (R),
    # which is the fatty acyl chain. Locate the neighbor of the carbonyl atom that is not the oxygen (or the sulfur).
    acyl_start_idx = None
    for neighbor in carbonyl_atom.GetNeighbors():
        # Exclude the oxygen (double-bonded to carbon).
        if neighbor.GetAtomicNum() == 8:
            continue
        # Exclude the sulfur atom from the thioester bond.
        if neighbor.GetIdx() == sulfur_idx:
            continue
        acyl_start_idx = neighbor.GetIdx()
        break
    if acyl_start_idx is None:
        return False, "No fatty acyl chain attached to the thioester carbonyl found"

    # Now we want to “collect” the fatty acyl chain atoms. A simple BFS starting from the acyl chain start,
    # while deliberately not going into the carbonyl atom (and so, not entering the CoA moiety), should serve.
    acyl_chain_atom_idxs = set()
    queue = [acyl_start_idx]
    while queue:
        current_idx = queue.pop(0)
        if current_idx in acyl_chain_atom_idxs:
            continue
        acyl_chain_atom_idxs.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            # Do not cross the bond back to the carbonyl (which connects to the CoA portion).
            if nbr.GetIdx() == carbonyl_idx:
                continue
            if nbr.GetIdx() not in acyl_chain_atom_idxs:
                queue.append(nbr.GetIdx())

    # With the set of atom indices assumed to belong to the fatty acyl chain, count the carbon-carbon double bonds.
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Consider only bonds where both atoms belong to the acyl chain.
        if a1.GetIdx() in acyl_chain_atom_idxs and a2.GetIdx() in acyl_chain_atom_idxs:
            # Check if the bond is a double bond and that both atoms are carbon (atomic number 6).
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                    double_bond_count += 1

    # According to the definition, the fatty acyl chain must contain exactly one C=C bond.
    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} C=C bond(s) instead of exactly one"

    # Optionally one might further try to validate that the remainder of the molecule indeed has features of CoA.
    # For simplicity we assume that if a thioester exists and the acyl chain qualifies, then the molecule is a fatty acyl-CoA.
    
    return True, "Molecule is a monounsaturated fatty acyl-CoA with one C=C bond in the fatty acyl chain"


# For testing purposes, you might run:
if __name__ == '__main__':
    test_smiles = "CCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    valid, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(valid, "->", reason)
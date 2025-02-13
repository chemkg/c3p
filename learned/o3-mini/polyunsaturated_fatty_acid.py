"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Polyunsaturated fatty acid (PUFA)
Definition: Any fatty acid containing more than one C=C double bond.
Acids in this group are reported to have cardioprotective effects;
and levels are lowered in chronic fatigue syndrome.

This implementation first requires that the molecule is a single fragment
(without salts) and that it contains at least one free (non‐esterified) carboxyl group.
We then count the number of carbon–carbon double bonds (ignoring carbonyl C=O).
Finally, we locate the “acyl chain” by following the neighbor of the carboxyl carbon
(which should be a carbon) and computing the longest contiguous carbon chain length.
We also require that the molecule has no rings.
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid (PUFA) based on its SMILES string.
    Criteria:
      (1) Molecule parses and is a single fragment (to avoid salts or complexes)
      (2) Contains at least one free (non‐esterified) carboxyl group (either protonated or deprotonated)
          – we do not insist on a strict “terminal” pattern.
      (3) Contains at least 2 carbon–carbon double bonds (ignoring C=O bonds in the acid group)
      (4) Contains no rings (a typical fatty acid is acyclic)
      (5) Has an acyl chain (the chain beginning as the neighbor of the acid carbon) of sufficient length.
         (We compute the longest contiguous carbon path from that neighbor.)
      (6) The elemental composition should be “fatty‐acid–like”, i.e. no extra atoms beyond C, O, H (except for a very few cases)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a PUFA, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Exclude molecules that come as salts/multiple fragments.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) != 1:
        return False, "Molecule contains multiple fragments (possible salt or conjugate)"
    
    # Optional: Check that the molecule is mostly C, H, and O.
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (1, 6, 8):
            return False, f"Atom {atom.GetSymbol()} found; not a typical fatty acid (only C, H, O are allowed)"
    
    # Reject molecules with rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a typical free fatty acid"
    
    # Look for free carboxyl groups (the acid functionality, not part of an ester).
    # We look for both protonated ([CX3](=O)[OX1H]) and deprotonated ([CX3](=O)[O-]) forms.
    carboxyl_prot = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    carboxyl_deprot = Chem.MolFromSmarts("[CX3](=O)[O-]")
    matches = mol.GetSubstructMatches(carboxyl_prot) + mol.GetSubstructMatches(carboxyl_deprot)
    if not matches:
        return False, "No free carboxyl group found"
        
    # For PUFA we assume only one carboxyl group is present.
    if len(matches) > 1:
        return False, "Multiple free carboxyl groups found; not a typical free fatty acid"
    
    # Identify the carboxyl carbon from the match.
    # In our SMARTS the first atom is the carboxyl carbon.
    carboxyl_match = matches[0]
    acid_carbon = mol.GetAtomWithIdx(carboxyl_match[0])
    
    # Now, relax the “terminal” requirement: we simply look for a carbon neighbor that is not oxygen.
    chain_start = None
    for nbr in acid_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            chain_start = nbr
            break
    if chain_start is None:
        return False, "Carboxyl carbon does not connect to an alkyl chain"
        
    # Count the number of carbon–carbon double bonds (ignore non-C=C bonds, such as acid C=O).
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom(); a2 = bond.GetEndAtom()
            # Count only C=C bonds.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1
    if double_bond_count < 2:
        return False, f"Only {double_bond_count} carbon–carbon double bond(s) found; need at least 2"
        
    # Define a helper function to compute the longest simple path (in terms of number of C atoms)
    # starting from a given carbon atom. Since PUFA chains are typically not very large,
    # a recursive DFS is acceptable.
    def longest_chain(atom, visited):
        max_len = 1  # count current atom
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                visited.add(nbr.GetIdx())
                length = 1 + longest_chain(nbr, visited)
                visited.remove(nbr.GetIdx())
                if length > max_len:
                    max_len = length
        return max_len

    # Calculate the acyl chain length starting from the chain_start atom.
    # We do not count the acid carbon.
    chain_length = longest_chain(chain_start, {chain_start.GetIdx()})
    if chain_length < 4:
        return False, f"Acyl chain too short (only {chain_length} carbons found); not a fatty acid"
    
    # (Optional) We could require that most carbon atoms in the molecule belong to the chain.
    # Count total carbon atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # The chain (plus acid carbon) should ideally be a high fraction of the carbons.
    if (chain_length + 1) / total_carbons < 0.6:
        return False, f"Acyl chain (with acid carbon) only accounts for {(chain_length+1)}/{total_carbons} carbons; unusual fatty acid structure"
        
    return True, (f"Contains a free carboxyl group, {double_bond_count} C=C bonds, and an acyl chain of "
                  f"{chain_length+1} carbons (including acid carbon), meets PUFA criteria")

# (Optional) simple tests; you can remove or modify these when integrating into a larger package.
if __name__ == "__main__":
    test_smiles = [
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",  # farnesoic acid (should be True)
        "OC(=O)C\\C=C\\C=C/C=C=CC#CC#C",  # mycomycin (should be True)
        "C(=C/CCCCC)\\C=C\\C(O)=O",       # (2E,4E)-deca-2,4-dienoic acid (should be True)
        "OC(CCCCCCCC(O)=O)C=CC(O)C(O)CC=CCC",  # 9,12,13-Trihydroxyoctadeca-10,15-dienoic acid (should be True)
        "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCC(O)=O",  # long chain PUFA (should be True)
    ]
    for smi in test_smiles:
        result, reason = is_polyunsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
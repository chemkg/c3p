"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: Straight-chain saturated fatty acid (defined as any saturated fatty acid lacking a side‐chain)

A valid straight-chain free fatty acid must have exactly one free carboxylic acid group 
(protonated or deprotonated) and an unbranched (linear) carbon chain. Only C, O, H (and D) are allowed.
Modifications such as hydroxyl or oxo groups are allowed if they do not constitute a side chain.
Molecules with additional acid groups, extra carbons (branches) or extra atoms (N, P, S, etc.) are rejected.
"""

from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid (i.e. has exactly one free carboxylic acid group 
    and an unbranched aliphatic carbon chain) based solely on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the definition.
        str: Explanation for the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude multi‐fragment species (e.g. salts)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (e.g. a salt)"
    
    # Allow only C, O, H and D (atomic numbers: 6, 8, 1, 2)
    allowed_atomic_nums = {1, 2, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains disallowed atom: {atom.GetSymbol()}"
    
    # Find free carboxylic acid group(s); allow both protonated (-C(=O)OH) and deprotonated (-C(=O)[O-]) forms.
    acid_smarts_list = ["C(=O)[OH]", "C(=O)[O-]"]
    acid_matches = []
    for smarts in acid_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(patt)
        # Append the matches (tuple of atom indices) if found.
        acid_matches.extend(matches)
    # Remove duplicates (if any)
    acid_matches = list(set(acid_matches))
    if len(acid_matches) != 1:
        return False, f"Expected exactly one free carboxylic acid group; found {len(acid_matches)}"
    
    # Choose the free acid match. The first atom in the match is the carboxyl carbon.
    acid_C_idx = acid_matches[0][0]
    
    # Build the main carbon chain: get all carbons reachable from the acid carbon via single, non-aromatic C–C bonds.
    chain_set = set()
    stack = [acid_C_idx]
    while stack:
        current = stack.pop()
        if current in chain_set:
            continue
        chain_set.add(current)
        atom = mol.GetAtomWithIdx(current)
        for bond in atom.GetBonds():
            # Only proceed along single, non-aromatic bonds
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.GetIsAromatic():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6:  # only follow carbon
                    if nbr.GetIdx() not in chain_set:
                        stack.append(nbr.GetIdx())
    
    # For a fatty acid, we require a minimum chain length.
    if len(chain_set) < 4:
        return False, "Fatty acid chain too short (less than 4 carbons)"
    
    # Verify the carbon subgraph forms a linear (unbranched) chain.
    # For each carbon in chain_set, count neighbors (within chain_set) connected by single bonds.
    adjacency = {idx: set() for idx in chain_set}
    for idx in chain_set:
        atom = mol.GetAtomWithIdx(idx)
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.GetIsAromatic():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in chain_set:
                    adjacency[idx].add(nbr.GetIdx())
    # In a linear chain, there should be exactly 2 endpoints with degree 1 and all internal carbons degree 2.
    endpoints = [idx for idx, neigh in adjacency.items() if len(neigh) == 1]
    internals = [idx for idx, neigh in adjacency.items() if len(neigh) == 2]
    extras = [idx for idx, neigh in adjacency.items() if len(neigh) not in (1, 2)]
    if not (len(endpoints) == 2 and len(internals) == (len(chain_set) - 2) and len(extras) == 0):
        return False, "Carbon chain is branched or not linear"
    
    # In addition, ensure that no carbon outside the chain is attached to any carbon in the chain.
    # This catches cases where extra alkyl side chains are appended to the main chain.
    for idx in chain_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in chain_set:
                return False, f"Branching detected: extra carbon substituent attached to carbon index {idx}"
    
    return True, "Molecule is a straight-chain saturated fatty acid lacking side-chains"

# Example usage (testing the function on selected SMILES):
if __name__ == '__main__':
    # Some examples from the outcomes:
    test_examples = [
        # True positives:
        ("CCCCCCCCCCCCCCCCCCCC(O)=O", "icosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "hexacosanoic acid"),
        ("C[C@@H](O)CCCCCCCCCCCCCCCCCCC(O)=O", "(20R)-20-hydroxyhenicosanoic acid"),
        ("OC(C)CCCCCCCCCCCCCCCCCCC(=O)O", "20-hydroxyhenicosanoic acid"),
        ("CCCC(O)=O", "butyric acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "heptacosanoic acid"),
        ("OCCCCCCCCCCCCCCCCCC(O)=O", "20-hydroxyicosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "triacontanoic acid"),
        ("CCCCCCCCCCCCCCCCCCC(O)=O", "nonadecanoic acid"),
        ("CCCCCCCC(O)=O", "octanoic acid"),
        # False positive examples:
        ("COC(=O)CCC(O)=O", "monomethyl succinate"),
        ("OC(CC([O-])=O)C([O-])=O", "malate(2-)"),
        ("OCCCCCCCCCCCCCCCCCCCC(O)=O", "20-hydroxyicosanoic acid (duplicate sample)"),
    ]
    for smi, name in test_examples:
        res, reason = is_straight_chain_saturated_fatty_acid(smi)
        outcome = "CORRECT" if res else "WRONGLY REJECTED"
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {res}\nReason: {reason}\n")
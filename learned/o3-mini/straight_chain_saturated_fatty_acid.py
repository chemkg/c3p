"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: Straight-chain saturated fatty acid

Definition:
    A straight-chain saturated fatty acid (free acid) is defined as a molecule that:
      • is a single fragment containing only the atoms C, O, H and D,
      • contains exactly one free carboxylic acid group (in either protonated –C(=O)OH or deprotonated –C(=O)[O-] form),
      • does not have any additional ester (or other acyl) groups,
      • has no unsaturation (i.e. the carbon–carbon bonds are single),
      • and all carbon atoms form one connected “skeleton” that is a simple path (that is, linear with exactly two ends)
             with the carboxyl carbon at one end.
Any extra carbon substituents (side‐chains) will make the carbon skeleton branched.
Modifications like hydroxyl or oxo groups are allowed provided they do not add extra carbons.
"""

from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid lacking side-chains,
    according to the definition given above.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise.
        str: Explanation for the result.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude multi‐fragment species (e.g. salts)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (e.g. a salt)"
    
    # Allow only C, O, H and D (atomic numbers 6, 8, 1, 2)
    allowed_atomic_nums = {1, 2, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains disallowed atom: {atom.GetSymbol()}"
    
    # Disallow ester (acyl) groups.
    # A free fatty acid must not include an ester (i.e. a pattern C(=O)O[C])
    ester_patt = Chem.MolFromSmarts("C(=O)O[C]")
    if mol.HasSubstructMatch(ester_patt):
        return False, "Molecule contains an ester group substituent"
    
    # Look for free carboxylic acid groups.
    # Accept both protonated (-C(=O)OH) and deprotonated (-C(=O)[O-]) forms.
    acid_smarts = ["C(=O)[OH]", "C(=O)[O-]"]
    acid_matches = []
    for patt in acid_smarts:
        submol = Chem.MolFromSmarts(patt)
        acid_matches.extend(mol.GetSubstructMatches(submol))
    # Remove duplicate matches (if any)
    acid_matches = list(set(acid_matches))
    if len(acid_matches) != 1:
        return False, f"Expected exactly one free carboxylic acid group; found {len(acid_matches)}"
    
    # The free acid group match: assume the first atom in the match is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    
    # Build the carbon skeleton: the subgraph of all carbon atoms connected by single, non‐aromatic bonds.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return False, "No carbon atoms found in molecule"
    
    # Build an undirected graph among carbon atoms:
    # graph: dict where key = carbon atom idx, value = set of neighboring carbon idx (only if connected by a SINGLE bond)
    carbon_graph = {idx: set() for idx in carbon_idxs}
    for bond in mol.GetBonds():
        # Only consider bonds between carbon atoms that are single bonds.
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                i1 = a1.GetIdx()
                i2 = a2.GetIdx()
                if i1 in carbon_graph and i2 in carbon_graph:
                    carbon_graph[i1].add(i2)
                    carbon_graph[i2].add(i1)
        else:
            # Disallow any unsaturation (i.e. non-single bonds) between carbons.
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Carbon skeleton shows unsaturation (non-single C–C bond)"
    
    # Check that the carbon graph is connected.
    visited = set()
    stack = [carbon_idxs[0]]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        for nbr in carbon_graph[current]:
            if nbr not in visited:
                stack.append(nbr)
    if set(carbon_idxs) != visited:
        return False, "Carbon skeleton is not fully connected, indicating fragmentation or extra substituents"
    
    # For a linear (straight-chain) carbon skeleton, the graph must be a simple path.
    # In a simple path (with n > 2), exactly two carbons (the terminal ones) have degree 1 and all others degree 2.
    endpoints = [idx for idx, neigh in carbon_graph.items() if len(neigh) == 1]
    internals = [idx for idx, neigh in carbon_graph.items() if len(neigh) == 2]
    extras = [idx for idx, neigh in carbon_graph.items() if len(neigh) not in (1, 2)]
    if len(carbon_idxs) < 2:
        return False, "Not enough carbons for a fatty acid chain"
    if not (len(endpoints) == 2 and len(internals) == (len(carbon_idxs) - 2) and len(extras) == 0):
        return False, "Carbon skeleton is branched (contains extra carbon substituents) or is not a simple linear chain"
    
    # In a free fatty acid the carboxyl carbon must be one of the terminal carbons.
    if acid_carbon_idx not in endpoints:
        return False, "The free acid carboxyl carbon is not at the end of the carbon chain"
    
    return True, "Molecule is a straight-chain saturated fatty acid lacking side-chains"

# Example usage (testing the function on selected SMILES):
if __name__ == '__main__':
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
        # True hydroxy fatty acids:
        ("OCCCCCCCCCCCCC(O)=O", "13-hydroxytridecanoic acid"),
        ("OCCCCCCCCCCCCCCC(O)=O", "15-hydroxypentadecanoic acid"),
        # Examples that should be rejected:
        ("COC(=O)CCC(O)=O", "monomethyl succinate"),
        ("OCC(=O)[C@@H](O)[C@H](O)[C@H](O)C([O-])=O", "keto-D-fructuronate"),
        ("[C@H]1(CCCCCCCC(=O)[O-])[C@@H](CCCCCCCCO)O1", "(9S,10R)-9,10-epoxy-18-hydroxyoctadecanoate"),
        ("CC(O)CCCCCCCCCCCCCC(O)=O", "15-hydroxypalmitic acid"),
        ("C1[C@]2([C@]3([C@@](C4=C(CC3)C=C(C=C4)O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C(O)=O)(CC[C@@]2(C(=O)[C@@H]1O)C)[H])[H])[H]", "16alpha-hydroxyestrone 3-O-beta-D-glucuronide"),
    ]
    for smi, name in test_examples:
        res, reason = is_straight_chain_saturated_fatty_acid(smi)
        outcome = "CORRECT" if res else "WRONGLY CLASSIFIED"
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {res}\nReason: {reason}\n")
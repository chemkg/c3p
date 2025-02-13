"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: CHEBI: Straight-chain saturated fatty acid (defined as any saturated fatty acid lacking a side‐chain)

A straight-chain fatty acid must have a single free carboxylic acid group and an unbranched carbon skeleton 
(with only permitted modifications such as hydroxyl or oxo groups attached to the chain).
Additionally, the molecule must consist basically only of carbon and oxygen (aside from hydrogens and possibly deuterium).
Molecules with extra heteroatoms (e.g. N, P, S), salts, or branches in the carbon chain are rejected.
"""

from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid (i.e. has exactly one free carboxylic acid group 
    and an unbranched carbon backbone) based solely on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a straight-chain saturated fatty acid, False otherwise.
        str: Explanation for the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude multi‐fragment species (e.g. salts)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (e.g. a salt)"
    
    # Require that the molecule contains only allowed atoms.
    # Here we allow carbon (6), oxygen (8), hydrogen (1) and deuterium (2). 
    # (Any other element such as N, P, S will cause a rejection.)
    allowed_atomic_nums = {1, 2, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains disallowed atom: {atom.GetSymbol()}"
    
    # Look for a carboxylic acid substructure.
    # Try both protonated acid (-C(=O)OH) and deprotonated (-C(=O)[O-]).
    acid_smarts1 = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_smarts1)
    if not acid_matches:
        acid_smarts2 = Chem.MolFromSmarts("C(=O)[O-]")
        acid_matches = mol.GetSubstructMatches(acid_smarts2)
        if not acid_matches:
            return False, "No carboxylic acid group found"
    # Pick the first match. The first atom in the match is the carboxyl carbon.
    acid_C_idx = acid_matches[0][0]
    
    # Build a carbon-only subgraph (indices of atoms that are carbon) which are connected by single bonds.
    # We also require these bonds not to be aromatic.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # We perform a DFS starting from the acid carbon.
    component = set()
    stack = [acid_C_idx]
    while stack:
        current = stack.pop()
        if current in component:
            continue
        component.add(current)
        atom = mol.GetAtomWithIdx(current)
        for bond in atom.GetBonds():
            # Only consider single bonds that are not aromatic.
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.GetIsAromatic():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in carbon_idxs:
                    if nbr.GetIdx() not in component:
                        stack.append(nbr.GetIdx())
    
    if acid_C_idx not in component:
        return False, "Carboxyl carbon not found in carbon skeleton"
    
    # We require a minimum chain length. (But note that very short acids like formic acid are not typical fatty acids.)
    # Here we demand at least 4 carbons (e.g. butyric acid = C4).
    if len(component) < 4:
        return False, "Fatty acid chain too short"
    
    # Next, check that the carbon skeleton (the connected component containing the acid carbon)
    # is a simple, unbranched chain. That means that when we count connections (edges) within this set,
    # there should be exactly 2 atoms with only one carbon neighbor (endpoints)
    # and every other carbon should have exactly 2 carbon neighbors.
    # Build an adjacency dictionary for carbons in the component.
    adjacency = {idx: set() for idx in component}
    for idx in component:
        atom = mol.GetAtomWithIdx(idx)
        for bond in atom.GetBonds():
            # Again, only consider single, non-aromatic bonds.
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.GetIsAromatic():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in component:
                    adjacency[idx].add(nbr.GetIdx())
    
    # Count degrees (number of carbon neighbors in the component)
    degrees = {idx: len(neighbs) for idx, neighbs in adjacency.items()}
    endpoints = [idx for idx, deg in degrees.items() if deg == 1]
    internals = [idx for idx, deg in degrees.items() if deg == 2]
    extras = [idx for idx, deg in degrees.items() if deg not in (1, 2)]
    
    if len(component) == 1:
        return False, "Only one carbon in skeleton"
    if not (len(endpoints) == 2 and len(internals) == (len(component) - 2) and len(extras) == 0):
        return False, "Branching or non-linear (e.g. cyclic) carbon chain detected"
    
    # If we reached here, we have exactly one unbranched carbon chain (including the acid carbon)
    return True, "Molecule is a straight-chain saturated fatty acid lacking side-chains"


# Example usage (testing the function on selected SMILES):
if __name__ == '__main__':
    test_examples = [
        # True positives (free acids with unbranched chains)
        ("CCCCCCCCCCCCCCCCCCCC(O)=O", "icosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "hexacosanoic acid"),
        ("C[C@@H](O)CCCCCCCCCCCCCCCCCCC(O)=O", "(20R)-20-hydroxyhenicosanoic acid"),
        ("OC(C)CCCCCCCCCCCCCCCCCCC(=O)O", "20-hydroxyhenicosanoic acid"),
        ("CCCC(O)=O", "butyric acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "heptacosanoic acid"),
        # False positives (should fail)
        ("O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CCC(=O)N", "Asp-Gln-Pro"),
        ("[Na+].CCCC([O-])=O", "sodium butyrate"),
    ]
    for smi, name in test_examples:
        res, reason = is_straight_chain_saturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {res}\nReason: {reason}\n")
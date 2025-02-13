"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: short-chain fatty acid
Definition:
  An aliphatic monocarboxylic acid with a chain length (including the carboxyl carbon) of less than C6.
  In addition, substituents on the backbone must be hydrocarbon in nature – if any non‐hydrocarbon substituent is present
  (apart from the carboxyl group itself), the compound is not normally regarded as a short‐chain fatty acid.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    
    A short-chain fatty acid in this context is defined as an aliphatic monocarboxylic acid
    (i.e. containing exactly one carboxyl group) whose longest continuous carbon chain (including
    the carboxyl carbon) has fewer than 6 carbon atoms. Furthermore, any substituent attached
    to the backbone that is not purely hydrocarbon (only carbon atoms, with implicit hydrogens) 
    disqualifies the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Molecule must be acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not purely aliphatic"
    
    # Identify the carboxyl group.
    # We use a SMARTS pattern for a carboxylic acid: a trigonal carbon with a double-bonded O and an -OH.
    acid_query = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_query)
    if len(acid_matches) != 1:
        return False, "Molecule does not have exactly one carboxyl group"
    
    # In the match, the first atom is the carboxyl carbon and the remaining atoms are the oxygens.
    acid_match = acid_matches[0]
    acid_carbon = acid_match[0]
    acid_oxygens = set(acid_match[1:])
    
    # Build a graph of all carbon atoms (by their indices) that are connected via C–C bonds.
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_graph = {idx: [] for idx in carbon_atoms}
    for idx in carbon_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                nidx = nbr.GetIdx()
                if nidx in carbon_graph:
                    carbon_graph[idx].append(nidx)
    
    # Restrict to the connected component of carbons reachable from the acid carbon.
    def dfs_collect(start):
        comp = set()
        stack = [start]
        while stack:
            current = stack.pop()
            if current not in comp:
                comp.add(current)
                for neigh in carbon_graph[current]:
                    if neigh not in comp:
                        stack.append(neigh)
        return comp

    acyl_component = dfs_collect(acid_carbon)
    if acid_carbon not in acyl_component:
        return False, "Acid carbon is not in the carbon component"
    
    # Next, find the longest simple path (in terms of number of carbon atoms) starting at the acid carbon.
    max_path = []
    def dfs_path(current, path, visited):
        nonlocal max_path
        if len(path) > len(max_path):
            max_path = path.copy()
        for neigh in carbon_graph[current]:
            if neigh not in visited:
                visited.add(neigh)
                path.append(neigh)
                dfs_path(neigh, path, visited)
                path.pop()
                visited.remove(neigh)
    dfs_path(acid_carbon, [acid_carbon], set([acid_carbon]))
    chain_length = len(max_path)
    
    # Check that the length (number of carbons in the backbone) is less than 6.
    if chain_length >= 6:
        return False, f"Longest carbon chain has {chain_length} carbons, exceeding the short-chain limit (<6)"
    
    # The chosen longest chain is our "backbone". (If there are more than one longest chain,
    # we choose one arbitrarily; note that if a heteroatom is attached to any backbone carbon
    # in any longest chain, it should disqualify the molecule.)
    backbone = set(max_path)

    # For each carbon atom in the backbone, examine substituents (neighbors not in the backbone).
    # For the acid carbon, ignore the oxygen atoms that are part of the acid group.
    # For each branch (which may be more than one atom long) we recursively check that the entire
    # branch is hydrocarbon. (Allowed atoms are carbon (6) or hydrogen (1); note that hydrogens are
    # usually implicit, so only explicit atoms need be checked.)
    def branch_is_hydrocarbon(start, visited):
        """DFS on the branch to ensure all atoms are carbon (if not hydrogen)"""
        if start in visited:
            return True
        visited.add(start)
        atom = mol.GetAtomWithIdx(start)
        # Only carbon is allowed in branches (if not simply hydrogen, which would be implicit)
        if atom.GetAtomicNum() != 6:
            return False
        for nbr in atom.GetNeighbors():
            # Do not traverse back into any atom that is in the backbone.
            if nbr.GetIdx() in backbone:
                continue
            if nbr.GetIdx() in visited:
                continue
            # Allow hydrogen if explicit (though typically they are implicit)
            if nbr.GetAtomicNum() not in (6, 1):
                return False
            if nbr.GetAtomicNum() == 6:
                if not branch_is_hydrocarbon(nbr.GetIdx(), visited):
                    return False
        return True

    for idx in backbone:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nidx = nbr.GetIdx()
            # For the acid carbon, skip the oxygens that are part of the acid.
            if idx == acid_carbon and nidx in acid_oxygens:
                continue
            # Skip if the neighbor is part of the backbone.
            if nidx in backbone:
                continue
            # If the neighbor is not carbon (or hydrogen), then it is a non-hydrocarbon substituent.
            if nbr.GetAtomicNum() not in (6, 1):
                return False, "Non-hydrocarbon substituent found attached to the acyl chain"
            # If the neighbor is carbon, traverse its branch.
            if nbr.GetAtomicNum() == 6:
                if not branch_is_hydrocarbon(nidx, set()):
                    return False, "Non-hydrocarbon substituent found attached to the acyl chain"
    return True, "Molecule is an aliphatic monocarboxylic acid with a carbon chain length of less than C6"

# (Optional) Main block to test a few SMILES strings.
if __name__ == "__main__":
    test_smiles = [
        "CCCCC(O)=O",                 # valeric acid (should be True)
        "CC(C)=CC(O)=O",              # 3-methylbut-2-enoic acid (should be True)
        "CCC(C)C(O)=O",               # (R)-2-methylbutyric acid (should be True)
        "CCCC(O)=O",                  # butyric acid (should be True)
        "CCCC(CC=C)C(O)=O",           # 2-n-Propyl-4-pentenoic acid (should be False: chain too long)
        "OC(=O)C(CC)(CC)C",           # 2-ethyl-2-methyl-butanoic acid (should be False: extra substituents)
        "CCC(CO)C(O)=O",              # 2-ethylhydracrylic acid (should be False: non-hydrocarbon substituent)
    ]
    for smi in test_smiles:
        result, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\n  Is short-chain fatty acid? {result}\n  Reason: {reason}\n")
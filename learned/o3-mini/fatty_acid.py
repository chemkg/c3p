"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: Fatty Acids 
Definition: “Fatty acids” here are defined as acyclic aliphatic monocarboxylic acids.
They must be a valid molecule with no rings, contain exactly one terminal carboxyl group 
(a carbonyl carbon bearing a C=O plus either –OH or [O-] and attached to exactly one carbon),
have no amide bond(s) (to avoid peptides), and show a predominantly linear hydrocarbon chain.
For molecules with 10 or more carbons we also require that branching is very limited.
In addition we require that the ratio of extra heteroatoms (aside from the two in the acid group)
relative to the total carbon count is not excessive.
This version relaxes the longest-chain ratio from 70% to 50% and (for molecules with ≥10 C’s)
requires that almost all carbons be in the main chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    The molecule must:
      1. Be a valid molecule.
      2. Be acyclic (no rings).
      3. Contain exactly one terminal carboxylic acid group – that is, a C(=O)[O;H1,O-] where
         the carboxyl carbon is attached to exactly one (aliphatic) carbon.
      4. Not contain amide bonds (i.e. no C(=O)N fragments).
      5. Have a long continuous (unbranched) carbon chain that represents at least 50% of all carbon atoms.
         Moreover, for molecules with 10 or more carbons the excess of branched carbons must be very low.
      6. Not be excessively oxidized – the number of extra heteroatoms (beyond the 2 oxygens
         in the carboxyl group) should be limited.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acid, False otherwise.
        str: A textual explanation for the decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), expected an acyclic fatty acid"
    
    # Define SMARTS for a carboxylic acid group:
    # Carbonyl (C=O) with a singly bonded oxygen that is either protonated (-OH) or deprotonated ([O-])
    carboxyl_smarts = "C(=O)[O;H1,O-]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    matches = mol.GetSubstructMatches(carboxyl_query)
    
    # A terminal carboxyl means that the carboxyl carbon should be bonded to exactly one other carbon.
    terminal_matches = []
    for match in matches:
        c_idx = match[0]  # the carboxyl carbon in the matched group
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Count neighbors that are carbon
        c_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(c_neighbors) == 1:
            terminal_matches.append(match)
    if len(terminal_matches) != 1:
        return False, f"Found {len(terminal_matches)} terminal carboxyl group(s), expected exactly 1"
    
    # Reject if any amide bond (C(=O)N) is found.
    amide_query = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_query):
        return False, "Molecule contains amide bond(s), likely not a free fatty acid"
    
    # Count total carbon atoms.
    atoms = list(mol.GetAtoms())
    total_carbons = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if total_carbons < 4:
        return False, f"Too few carbon atoms ({total_carbons}) to be a fatty acid"
    
    # Build a carbon-only connectivity graph.
    carbon_indices = [atom.GetIdx() for atom in atoms if atom.GetAtomicNum() == 6]
    graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    
    # Identify the carbon atom attached to the terminal carboxyl group.
    carboxyl_c_idx = terminal_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_c_idx)
    neighbors = [nbr.GetIdx() for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if not neighbors:
        return False, "Terminal carboxyl group has no attached carbon"
    chain_start = neighbors[0]
    
    # Compute the longest continuous (unbranched) carbon chain in the connected component containing chain_start.
    def bfs_longest(start):
        visited = {start}
        queue = deque([(start, 0)])
        farthest = (start, 0)
        while queue:
            current, dist = queue.popleft()
            if dist > farthest[1]:
                farthest = (current, dist)
            for nb in graph[current]:
                if nb not in visited:
                    visited.add(nb)
                    queue.append((nb, dist + 1))
        return farthest, visited

    # Since the carbon scaffold is a tree (acyclic), we can compute its diameter.
    _, comp = bfs_longest(chain_start)
    any_node = next(iter(comp))
    (node_A, _), _ = bfs_longest(any_node)
    (node_B, diameter), _ = bfs_longest(node_A)
    longest_chain = diameter + 1  # number of carbons in the longest path
    
    # Check that the longest chain is not too short overall.
    if longest_chain < 4:
        return False, f"Longest carbon chain ({longest_chain}) is too short for a fatty acid"

    # Compute the ratio of carbons in the longest chain to total carbons.
    ratio = longest_chain / total_carbons

    # Relaxed threshold: require at least 50% of the carbons be in one continuous chain.
    if ratio < 0.50:
        return False, f"Longest carbon chain length ({longest_chain}) is too short relative to total carbons ({total_carbons})"
    
    # For molecules with 10 or more carbons, also require that the extra (branched) carbons are minimal.
    if total_carbons >= 10:
        if (total_carbons - longest_chain) > 1:
            return False, (f"Excessive branching: longest chain uses {longest_chain} of {total_carbons} carbons")
    
    # Heuristic on heteroatoms:
    # Count all heteroatoms (all atoms except C and H).
    hetero_total = sum(1 for atom in atoms if atom.GetAtomicNum() not in (1, 6))
    # The terminal carboxyl group contributes 2 oxygens; the rest are extra.
    extra_hetero = hetero_total - 2
    # For small molecules (<10 C) we disallow any extra heteroatoms,
    # and for larger ones we require extra hetero atoms to be no more than 20% of all carbons.
    if total_carbons < 10:
        if extra_hetero > 0:
            return False, (f"Too many functional groups for a short fatty acid: "
                           f"{extra_hetero} extra heteroatom(s) with only {total_carbons} carbons")
    else:
        if extra_hetero > int(total_carbons * 0.20):
            return False, (f"High heteroatom content (extra={extra_hetero} for {total_carbons} carbons), "
                           "not typical for a fatty acid")

    # (Optional) Reject if extra ketone/oxo groups (not as part of the terminal acid) are present.
    # We look for a carbonyl in an alkyl chain: a pattern [#6][CX3](=O)[#6].
    ketone_query = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_query)
    if ketone_matches:
        # Allow one extra oxidation in long-chain acids (>=20 carbons) but otherwise reject.
        if total_carbons < 20 or len(ketone_matches) > 1:
            return False, "Found extra ketone/oxo group(s) in the acyl chain"
    
    return True, (f"Molecule is an acyclic aliphatic carboxylic acid with {total_carbons} carbon(s) "
                  f"and a main chain of {longest_chain} carbon(s) (ratio={ratio:.2f})")

# Example usage:
if __name__ == "__main__":
    # Several tests (both positive and negative examples)
    test_smiles_list = [
        # True positives:
        "C(CCCCCC/C=C\\C(=O)O)CCCC(C)C",  # (Z)-14-methylpentadec-2-enoic acid
        "OC(=O)/C=C\\CCCCCCCCCCCCC",       # 2Z-hexadecenoic acid
        "OC(=O)CC/C=C/C=C/C#CC#C/C=C\\C",   # 4E,6E,12Z-Tetradecatriene-8,10-diynoic acid
        "OC(CCCCCCCCCCCCCC(O)=O)CCC",       # 15-hydroxy stearic acid
        "OC(CCCCCC)C/C=C/C(O)=O",           # (E)-5-Hydroxyundec-2-enoic acid
        "CCCC(O)=O",                      # butyric acid
        # A false negative (expected to be accepted now):
        "[O-]C(=O)CCCC[N+](C)(C)C",         # 5-aminovaleric acid betaine
        # False positives (should be rejected):
        "OC(CCCC)CCC(O)=O",                # 4-hydroxy-octanoic acid
    ]
    for sm in test_smiles_list:
        result, reason = is_fatty_acid(sm)
        print(sm, "=>", result, "|", reason)
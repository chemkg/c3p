"""
Classifies: CHEBI:28966 chlorophyll
"""
#!/usr/bin/env python3
"""
Classifies: Chlorophyll – magnesium porphyrins with a fifth ring and a long phytol side chain.
The criteria used are:
  1. The molecule must contain magnesium.
  2. At least one magnesium must be bound to four nitrogen atoms (each in a ring) – a proxy for the porphyrin.
  3. The molecule must have at least five rings (rough proxy for the fused porphyrin plus the extra ring).
  4. The molecule must have a long (acyclic) alkyl chain (at least 8 carbons) representing the phytol chain.
If any of these aren’t met, a reason is given.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determine if a molecule belongs to the chlorophyll family based on its SMILES string.

    The used criteria are:
       - The molecule contains a magnesium center.
       - That magnesium is coordinated to four nitrogen atoms that are part of rings (indicative of a porphyrin core).
       - The molecule has at least 5 rings (porphyrin macrocycle plus an extra “fifth” ring).
       - It possesses a long, mostly unbranched alkyl chain (at least 8 C atoms) that represents the phytol side chain.

    Args:
        smiles (str): SMILES string to be tested.

    Returns:
        (bool, str): True and a success message if it matches the chlorophyll criteria, or
                     False and a reason why not.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for magnesium atoms.
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Mg"]
    if not mg_atoms:
        return False, "No magnesium present: not a Mg-porphyrin"

    # 2. For at least one Mg, check that it is bound to 4 nitrogen atoms that are part of rings.
    mg_in_porhyrin = False
    for mg in mg_atoms:
        n_neighbors = [nbr for nbr in mg.GetNeighbors() if nbr.GetSymbol() == "N" and nbr.IsInRing()]
        if len(n_neighbors) == 4:
            mg_in_porhyrin = True
            break
    if not mg_in_porhyrin:
        return False, "Magnesium center not bound to 4 ring-nitrogen atoms (porphyrin core missing)"

    # 3. Count ring systems. We require at least 5 rings. (The porphyrin system normally appears as several fused rings.)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 5:
        return False, f"Only {len(rings)} rings found; expected at least 5 (porphyrin plus extra ring)"
    
    # 4. Detect a long alkyl chain (phytol chain). We search for a chain of non-ring carbon atoms.
    # Build a graph of non-ring carbon atoms only.
    nodes = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and not atom.IsInRing():
            idx = atom.GetIdx()
            nodes[idx] = []
    
    # Build connectivity among these non-ring carbons.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        if idx1 in nodes and idx2 in nodes:
            nodes[idx1].append(idx2)
            nodes[idx2].append(idx1)
    
    # Find the longest simple path in this (acyclic) graph of non-ring carbons.
    def dfs(node, visited):
        best = 1
        for neighbor in nodes[node]:
            if neighbor not in visited:
                path_len = 1 + dfs(neighbor, visited | {neighbor})
                if path_len > best:
                    best = path_len
        return best
    
    longest_chain = 0
    for node in nodes:
        chain_len = dfs(node, {node})
        if chain_len > longest_chain:
            longest_chain = chain_len

    # If the longest chain is less than 8 carbon atoms, we assume the phytol chain is missing.
    if longest_chain < 8:
        return False, f"Longest acyclic carbon chain is only {longest_chain} atoms long, expected at least 8 for the phytol chain"
    
    # If all tests pass, we consider the molecule as chlorophyll.
    return True, "Molecule has a magnesium-centered porphyrin core with an additional (fifth) ring and a long phytol chain"

# If desired, one might test the function with an example SMILES:
if __name__ == "__main__":
    example_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"  # chlorophyll a
    result, reason = is_chlorophyll(example_smiles)
    print(result, reason)
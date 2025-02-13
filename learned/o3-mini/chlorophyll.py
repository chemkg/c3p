"""
Classifies: CHEBI:28966 chlorophyll
"""
#!/usr/bin/env python3
"""
Classifies: Chlorophyll – magnesium porphyrins with an extra (fifth) ring and a long phytol chain.
Criteria used are:
  1. The molecule must have a magnesium atom.
  2. At least one Mg should be coordinated to at least 3 nitrogen atoms that lie in rings 
     (a relaxed proxy for a porphyrin macrocycle).
  3. The molecule must contain at least 5 rings (the porphyrin system [nominally 4 rings] plus one extra ring).
  4. The molecule must have a long, acyclic alkyl chain (at least 8 connected non‐ring carbon atoms)
     that represents the phytol side chain.
  5. In order to avoid mis‐classifying ionic derivatives (e.g. “chlorophyll a(1-)”), we reject molecules 
     with net negative charge.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule belongs to the chlorophyll family based on its SMILES string.

    The criteria are:
        - Contains a magnesium atom.
        - At least one magnesium atom is coordinated to at least 3 nitrogen atoms that are part
          of rings (serving as a proxy for the porphyrin core).
        - The molecule has at least 5 rings (porphyrin core plus an extra “fifth” ring).
        - It possesses a long acyclic alkyl chain (a chain of at least 8 non‐ring carbons) representing
          the phytol side chain.
        - The overall molecule must not be an anion.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the chlorophyll criteria, False otherwise.
        str: Explanation of the classification outcome.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if the overall molecule has a net negative charge.
    net_charge = Chem.GetFormalCharge(mol)
    if net_charge < 0:
        return False, f"Molecule has net charge {net_charge}; likely a negatively charged derivative."

    # 1. Check for the presence of magnesium.
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Mg"]
    if not mg_atoms:
        return False, "No magnesium present: not a magnesium porphyrin"
    
    # 2. Check that at least one Mg is coordinated to at least 3 nitrogen atoms that are in rings.
    mg_in_porhyrin = False
    for mg in mg_atoms:
        # Get neighbors that are nitrogen and are in a ring.
        n_neighbors = [nbr for nbr in mg.GetNeighbors() if nbr.GetSymbol() == "N" and nbr.IsInRing()]
        if len(n_neighbors) >= 3:
            mg_in_porhyrin = True
            break
    if not mg_in_porhyrin:
        return False, "Magnesium center not bound to at least 3 ring-nitrogen atoms (porphyrin core missing)"
    
    # 3. Count the total number of rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 5:
        return False, f"Only {len(rings)} rings found; expected at least 5 (four for the porphyrin core plus one extra ring)"

    # 4. Determine the length of the longest acyclic carbon chain.
    # We consider only non-ring carbons.
    nonring_carbons = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and not atom.IsInRing():
            idx = atom.GetIdx()
            nonring_carbons[idx] = []
    
    # Build connectivity among these non-ring carbons.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        if idx1 in nonring_carbons and idx2 in nonring_carbons:
            nonring_carbons[idx1].append(idx2)
            nonring_carbons[idx2].append(idx1)
    
    # DFS function to find longest simple path in the graph of non-ring carbons.
    def dfs(node, visited):
        best = 1
        for neighbor in nonring_carbons[node]:
            if neighbor not in visited:
                path_len = 1 + dfs(neighbor, visited | {neighbor})
                if path_len > best:
                    best = path_len
        return best
    
    longest_chain = 0
    for node in nonring_carbons:
        chain_len = dfs(node, {node})
        if chain_len > longest_chain:
            longest_chain = chain_len

    if longest_chain < 8:
        return False, f"Longest acyclic carbon chain is only {longest_chain} atoms long, expected at least 8 for the phytol side chain"
    
    return True, "Molecule has a magnesium-centered porphyrin core (with at least 3 coordinated ring N atoms), at least 5 rings overall, and a long phytol chain."

# For testing, you may uncomment the code block below.
if __name__ == "__main__":
    # Example: chlorophyll a SMILES (one known true positive)
    example_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, reason = is_chlorophyll(example_smiles)
    print(result, reason)
"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: Oligosaccharide
Definition: A compound in which monosaccharide units are joined by glycosidic linkages.
A valid oligosaccharide here should have at least two sugar rings that are linked via an ether (glycosidic) bond.
Note: This is a heuristic approach and may not cover all edge cases.
"""

from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    The heuristic used is:
       1. Identify sugar rings: rings of size 5 (furanose) or 6 (pyranose) that contain exactly one oxygen atom.
       2. Check that at least two sugar rings are detected.
       3. Verify that the sugar rings are connected via a glycosidic bond. Here we define a glycosidic link
          as either a direct bond between an atom in one sugar ring and an atom in another, or a two-bond path
          using an exocyclic oxygen that is not part of either ring.

    Args:
         smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is an oligosaccharide, False otherwise
         str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information (list of tuples of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()

    # Heuristic: identify sugar rings as rings of size 5 or 6 that contain exactly one oxygen atom.
    sugar_rings = []
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_count += 1
        # In a typical sugar ring, for a 5-membered ring (furanose), expect 1 O and 4 C;
        # For a 6-membered ring (pyranose), expect 1 O and 5 C.
        if len(ring) == 5 and oxygen_count == 1 and carbon_count == 4:
            sugar_rings.append(set(ring))
        elif len(ring) == 6 and oxygen_count == 1 and carbon_count == 5:
            sugar_rings.append(set(ring))

    if len(sugar_rings) < 2:
        return False, f"Fewer than 2 sugar rings detected (found {len(sugar_rings)})."

    # Build a connectivity graph between sugar rings.
    # Each node represents a sugar ring (by index in sugar_rings).
    # We add an edge if two sugar rings are connected by a glycosidic bond.
    #
    # Heuristic for glycosidic bond:
    # 1. Direct connection: if there is a bond directly between any atom in ring A and any atom in ring B.
    # 2. Bridging oxygen: if an atom 'a' in ring A is bonded to an oxygen (that is not part of ring A)
    #    which in turn is bonded to an atom 'b' that is in ring B.
    n = len(sugar_rings)
    graph = {i: set() for i in range(n)}

    # Iterate over pairs of sugar rings
    for i in range(n):
        for j in range(i+1, n):
            connected = False
            # Check direct bonds: iterate over bonds in molecule
            for bond in mol.GetBonds():
                a_idx = bond.GetBeginAtomIdx()
                b_idx = bond.GetEndAtomIdx()
                if a_idx in sugar_rings[i] and b_idx in sugar_rings[j]:
                    connected = True
                    break
                if a_idx in sugar_rings[j] and b_idx in sugar_rings[i]:
                    connected = True
                    break
            # If no direct bond, try to find a bridging oxygen.
            if not connected:
                # For each atom in ring i, look for an exocyclic neighbor that is an oxygen,
                # then see if that oxygen is bonded to an atom in ring j.
                for a_idx in sugar_rings[i]:
                    atom_a = mol.GetAtomWithIdx(a_idx)
                    for nbr in atom_a.GetNeighbors():
                        # Skip if neighbor is within the same ring
                        if nbr.GetIdx() in sugar_rings[i]:
                            continue
                        if nbr.GetAtomicNum() != 8:
                            continue
                        # Check neighbors of the oxygen
                        for nbr2 in nbr.GetNeighbors():
                            # Do not go back to the same atom
                            if nbr2.GetIdx() == a_idx:
                                continue
                            if nbr2.GetIdx() in sugar_rings[j]:
                                connected = True
                                break
                        if connected:
                            break
                    if connected:
                        break
            if connected:
                graph[i].add(j)
                graph[j].add(i)

    # Check if the sugar ring connectivity graph is connected.
    # Do a simple DFS starting from node 0.
    visited = set()
    stack = [0]
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        for nbr in graph[node]:
            if nbr not in visited:
                stack.append(nbr)

    if len(visited) != n:
        return False, "Sugar rings detected but not all are connected by glycosidic linkages."
    
    return True, "Found multiple sugar rings connected via glycosidic linkages; likely an oligosaccharide."

# Example usage:
if __name__ == '__main__':
    # alpha-mannobiose (disaccharide) example:
    test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"
    result, reason = is_oligosaccharide(test_smiles)
    print(result, reason)
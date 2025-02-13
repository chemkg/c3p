"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid
with sulfuric acid.
A valid steroid sulfate must contain:
  1. A steroid nucleus (approximated by a fused ring system of at least 4 rings with a high carbon ratio).
  2. A sulfate ester group (–OS(=O)(=O)[O or O-]) attached via an oxygen atom to one of the steroid carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Steps:
      1. Parse the SMILES string to a molecule.
      2. Identify a steroid-like fused ring system: we look for a cluster of at least 4 fused rings,
         and we require that most of the atoms in these rings are carbons.
      3. Search for a sulfate ester group. Instead of using an overly strict SMARTS pattern,
         we iterate through oxygen atoms bonded to a sulfur atom (candidate sulfate group).
         For the candidate S atom, we check that it is attached to at least three oxygens (expected for a sulfate)
         and then check that the oxygen attached to S is also bonded to a carbon that is part of the steroid nucleus.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a steroid sulfate, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # ------ STEP 1: Identify the steroid nucleus ------
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples representing atom indices in each ring
    if not atom_rings:
        return False, "No rings found in the molecule; cannot be a steroid."
    
    # Build a graph where each node is a ring (by index in atom_rings).
    # Two rings are fused if they share at least two atoms.
    ring_graph = {i: set() for i in range(len(atom_rings))}
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            if len(set(atom_rings[i]).intersection(atom_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected clusters (components) of fused rings.
    visited = set()
    fused_clusters = []
    for i in range(len(atom_rings)):
        if i not in visited:
            stack = [i]
            cluster = set()
            while stack:
                curr = stack.pop()
                if curr not in visited:
                    visited.add(curr)
                    cluster.add(curr)
                    stack.extend(ring_graph[curr] - visited)
            fused_clusters.append(cluster)
    
    # Now, look for a cluster with at least 4 fused rings that has a high ratio of carbon atoms.
    steroid_core_atoms = None
    for cluster in fused_clusters:
        if len(cluster) >= 4:
            atoms_in_cluster = set()
            for ring_idx in cluster:
                atoms_in_cluster.update(atom_rings[ring_idx])
            # Count carbons among these atoms.
            num_carbons = sum(1 for idx in atoms_in_cluster if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            ratio_carbons = num_carbons / len(atoms_in_cluster)
            if ratio_carbons >= 0.7:
                steroid_core_atoms = atoms_in_cluster
                break
    
    if steroid_core_atoms is None:
        return False, "No steroid-like fused ring system (four fused rings with mostly carbons) found."
    
    # ------ STEP 2: Find a sulfate ester group attached to the steroid nucleus ------
    # We iterate over all oxygen atoms that are candidates for the ester oxygen.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:  # only consider oxygen atoms
            continue
        # Check neighbors of the oxygen: we want one neighbor to be sulfur (candidate sulfate group)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # sulfur
                sulfate = nbr
                # For a sulfate group, we expect at least three oxygen neighbors.
                oxy_neighbors = [n for n in sulfate.GetNeighbors() if n.GetAtomicNum() == 8]
                if len(oxy_neighbors) < 3:
                    continue  # not enough oxygens to be a sulfate
                # Now, check if the current oxygen (the one bridging to sulfate) 
                # is also bonded to a carbon from the steroid nucleus.
                for second_nbr in atom.GetNeighbors():
                    if second_nbr.GetAtomicNum() == 6 and second_nbr.GetIdx() in steroid_core_atoms:
                        # We have found an oxygen attached to sulfur (sulfate)
                        # and also attached to a steroid carbon.
                        return True, "Molecule contains a steroid nucleus with a sulfate ester attached via a hydroxy group."
    return False, "No sulfate ester group attached to the steroid nucleus found."


# For simple testing, you can run this script directly.
if __name__ == "__main__":
    # Example test SMILES: you may substitute one of the given steroid sulfate examples.
    test_smiles = "O(S(O)(=O)=O)[C@@H]1CC=2[C@]([C@]3(CC[C@]4([C@]([C@@]3(CC2)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])(C)CC1"
    result, reason = is_steroid_sulfate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
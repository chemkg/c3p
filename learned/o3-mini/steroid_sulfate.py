"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid
with sulfuric acid.

A valid steroid sulfate must contain:
  1. A steroid nucleus – approximated here as a fused-ring system of at least 4 rings,
     in which at least 60% of the atoms are carbons.
  2. A sulfate ester group: an oxygen atom that bridges a carbon to a sulfate S (S bonded to at least 3 oxygens).
     We further require that the bridging oxygen is only connected to exactly two atoms (one carbon and one sulfur)
     and that the carbon is either in the nucleus or directly adjacent to an atom in the nucleus.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.

    Steps:
      1. Parse the SMILES string.
      2. Identify a steroid nucleus based on fused rings:
         - Obtain all rings and build a connectivity graph for rings sharing at least 2 atoms.
         - Look for a connected set of at least 4 rings.
         - In that set, require that at least 60% of the atoms are carbons.
      3. Look for a sulfate ester group:
         - For each oxygen atom that has exactly 2 bonds,
           check that one neighbor is sulfur (atomic number 16) and the other is a carbon (atomic number 6).
         - Verify that the sulfur (candidate sulfate) has at least three oxygen neighbors.
         - Finally, check that the carbon linked via the oxygen is either directly in the steroid nucleus
           or at least one of its neighbors is in the nucleus.
    
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
    
    # Build a graph where nodes represent rings.
    # Two rings are fused if they share at least 2 atoms.
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
    
    steroid_core_atoms = None
    # Look for a cluster of at least 4 rings that has a high carbon ratio.
    for cluster in fused_clusters:
        if len(cluster) >= 4:
            atoms_in_cluster = set()
            for ring_idx in cluster:
                atoms_in_cluster.update(atom_rings[ring_idx])
            num_carbons = sum(1 for idx in atoms_in_cluster if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            ratio_carbons = num_carbons / float(len(atoms_in_cluster))
            if ratio_carbons >= 0.6:
                steroid_core_atoms = atoms_in_cluster
                break
    if steroid_core_atoms is None:
        return False, "No steroid-like fused ring system (>=4 rings with at least 60% carbons) found."
    
    # Helper: Check if a given carbon atom is considered part of or adjacent to steroid nucleus.
    def is_in_or_adjacent_to_nucleus(carbon_atom):
        idx = carbon_atom.GetIdx()
        if idx in steroid_core_atoms:
            return True
        for nbr in carbon_atom.GetNeighbors():
            # if one neighbor is in the nucleus, consider it attached.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in steroid_core_atoms:
                return True
        return False

    # ------ STEP 2: Find a valid sulfate ester group attached to the steroid ------
    # Look for candidate oxygen atoms that are "bridging" between sulfur and carbon.
    for atom in mol.GetAtoms():
        # Only consider oxygen atoms.
        if atom.GetAtomicNum() != 8:
            continue
        # We require that this oxygen has exactly 2 bonds (one to a sulfate S and one to a carbon)
        if atom.GetDegree() != 2:
            continue
        neighbors = list(atom.GetNeighbors())
        # Identify if one neighbor is sulfur and the other is carbon.
        sulfur_neighbor = None
        carbon_neighbor = None
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 16:
                sulfur_neighbor = nbr
            elif nbr.GetAtomicNum() == 6:
                carbon_neighbor = nbr
        if sulfur_neighbor is None or carbon_neighbor is None:
            continue

        # Verify that the sulfur indeed has at least three oxygen neighbors.
        oxy_neighbors = [n for n in sulfur_neighbor.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 3:
            continue

        # Finally check if the bridging carbon is (or is close to) part of the steroid nucleus.
        if is_in_or_adjacent_to_nucleus(carbon_neighbor):
            return True, ("Molecule contains a steroid nucleus with a sulfate ester attached via a hydroxy group "
                          "(bridging O with degree=2 linking to a sulfate S and steroid or adjacent C).")
    
    return False, "No sulfate ester group attached to the steroid nucleus found."


# For testing
if __name__ == "__main__":
    # Example test SMILES; you may substitute one of the steroid sulfate examples.
    test_smiles = "C[C@]12CCC3C(C1CCC2=O)CCC4=C3C=CC(=C4)OS(=O)(=O)O"
    result, reason = is_steroid_sulfate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
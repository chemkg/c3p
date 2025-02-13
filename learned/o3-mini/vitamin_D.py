"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: Vitamin D compounds (fat‐soluble hydroxy seco‐steroids)

Definition:
  “Any member of a group of fat‐soluble hydroxy seco‐steroids that exhibit biological activity
   against vitamin D deficiency. Vitamin D can be obtained from sun exposure, food and supplements
   and is biologically inactive and converted into the biologically active calcitriol via double
   hydroxylation in the body.”

Heuristic criteria in this version:
  1. The SMILES must parse to a valid molecule.
  2. The molecular weight should be between ~250 and ~600 Da.
  3. The molecule must contain at least one hydroxyl (-OH) group.
  4. The molecule must contain a fused ring cluster (with 2–5 rings) whose atoms are mostly non‐aromatic.
  5. There must be at least one extended aliphatic side chain (i.e. a path of at least 3 carbon atoms)
     that is attached to an atom in the fused ring cluster and that ends in a hydroxyl group.
     
This version improves the side chain search by “tracing back” from hydroxyl groups external to the cluster.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Heuristic criteria:
      - Molecule must be valid.
      - Molecular weight between ~250 and ~600 Da.
      - Contains at least one hydroxyl (-OH) group.
      - Contains a fused ring cluster (2–5 rings) that is largely non‐aromatic.
      - Has an extended aliphatic side chain (chain of >=3 carbon atoms) attached to the fused ring
        cluster that terminates in a hydroxyl group (outside the fused cluster).
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple of (True, reason) if classified as vitamin D,
                   or (False, reason) otherwise.
    """
    # Parse SMILES into molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (between ~250 and ~600 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.1f} Da is out of vitamin D range (250-600 Da)"
    
    # Check for at least one hydroxyl group (-OH).
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl (-OH) group found"
    
    # Obtain all rings.
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Build a graph of rings. Two rings are connected if they share >=2 atoms.
    n_rings = len(all_rings)
    adjacency = {i: set() for i in range(n_rings)}
    ring_sets = [set(ring) for ring in all_rings]
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Find connected components (clusters) of rings using DFS.
    visited = set()
    clusters = []
    for i in range(n_rings):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                curr = stack.pop()
                if curr in visited:
                    continue
                visited.add(curr)
                comp.add(curr)
                stack.extend(adjacency[curr])
            clusters.append(comp)
    
    # Select the largest fused ring cluster.
    largest_cluster = max(clusters, key=lambda comp: len(comp))
    cluster_ring_count = len(largest_cluster)
    if cluster_ring_count < 2 or cluster_ring_count > 5:
        return False, f"Fused ring cluster has {cluster_ring_count} rings; expected between 2 and 5 for a vitamin D secosteroid core"
    
    # Get all atom indices that are part of the fused ring cluster.
    cluster_atoms = set()
    for idx in largest_cluster:
        cluster_atoms.update(ring_sets[idx])
    
    # Check that the fused ring cluster is mostly non‐aromatic.
    arom_count = 0
    for atom_idx in cluster_atoms:
        if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
            arom_count += 1
    aromatic_ratio = arom_count / len(cluster_atoms)
    if aromatic_ratio > 0.5:
        return False, f"Fused ring cluster is predominantly aromatic (ratio {aromatic_ratio:.2f}); not consistent with vitamin D secosteroid core"
    
    # Now, search for an extended aliphatic side chain (with >= 3 carbons) attached to the fused cluster.
    # We look for hydroxyl groups (-OH) that are not part of the fused ring cluster.
    # For each such hydroxyl, we attempt to trace an aliphatic (carbon-only) path from a carbon attached to
    # the hydroxyl back to the fused cluster.
    def path_to_cluster_from(start_idx, visited_local):
        """
        Recursive DFS that traverses only through carbon atoms (atomic number 6) not in cluster_atoms.
        Returns the number of carbon atoms counted along the path if a path connecting to any fused cluster
        atom is found; otherwise returns None.
        """
        # Check neighbors for connection to fused cluster.
        atom = mol.GetAtomWithIdx(start_idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in cluster_atoms:
                return 1  # path length of one (current carbon will be counted in caller)
        # Otherwise, traverse neighbors that are carbons and not visited.
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in cluster_atoms:
                continue
            if nb.GetAtomicNum() != 6:
                continue
            if nb_idx in visited_local:
                continue
            visited_local.add(nb_idx)
            ret = path_to_cluster_from(nb_idx, visited_local)
            if ret is not None:
                return ret + 1
        return None

    branch_found = False
    # For every hydroxyl oxygen not in cluster, check if it is part of an extended side chain.
    hydroxyl_atoms = mol.GetSubstructMatches(hydroxyl_pattern)
    for match in hydroxyl_atoms:
        oh_idx = match[0]
        if oh_idx in cluster_atoms:
            continue  # skip hydroxyls that are part of the fused cluster
        oh_atom = mol.GetAtomWithIdx(oh_idx)
        # Check that this oxygen actually is -OH (attached to at least one hydrogen)
        has_hydrogen = any(nb.GetAtomicNum() == 1 for nb in oh_atom.GetNeighbors())
        if not has_hydrogen:
            continue
        # For each neighbor of the hydroxyl oxygen that is carbon and not in the cluster,
        # try to trace back to the fused cluster.
        for nb in oh_atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb.GetAtomicNum() != 6:
                continue
            if nb_idx in cluster_atoms:
                continue  # if the hydroxyl is directly attached to the cluster, it is not considered an extended chain.
            # Use DFS to find a path back to the cluster. Count the number of carbons in the path.
            visited_local = set([nb_idx])
            path_length = path_to_cluster_from(nb_idx, visited_local)
            if path_length is not None and path_length >= 3:
                branch_found = True
                break
        if branch_found:
            break

    if not branch_found:
        return False, "No extended aliphatic side chain (with at least 3 carbons and a hydroxyl) attached to the fused ring cluster"
    
    # If all checks pass, classify as vitamin D.
    return True, "Matches vitamin D criteria: fused non‐aromatic secosteroid core with hydroxyl and appropriate side chain"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with one example: 24-fluoro-25-hydroxyvitamin D3
    test_smiles = "FC(CC[C@H]([C@@H]1[C@@]2([C@@H](CC1)/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)C)C)C(O)(C)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)
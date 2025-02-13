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

Heuristic criteria used in this improved version:
  1. Valid molecule from SMILES.
  2. Molecular weight between ~250 and ~600 Da.
  3. Contains at least one hydroxyl (-OH) group.
  4. Contains a fused ring cluster (2–5 rings) whose atoms are mostly non‐aromatic.
  5. Has at least one extended aliphatic side chain (≥3 connected carbons) attached to a ring
     of the fused cluster that also bears a hydroxyl group.
  
Note:
  This heuristic is an approximation. Many vitamin D compounds have “seco‐steroid” cores
  (with a broken or “cut” ring) and pendant chains that carry –OH groups. Many non‐vitamin D 
  compounds may have a fused ring system, so we attempt to reduce false positives by requiring 
  the extra chain feature.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Heuristic criteria:
      - Valid molecule that can be parsed from the SMILES.
      - Molecular weight between ~250 Da and ~600 Da.
      - Contains at least one hydroxyl group.
      - Contains a fused ring cluster (2–5 rings) where most atoms are non‐aromatic.
      - Has an aliphatic side chain (≥~3 carbons) attached to the fused ring cluster that 
        contains at least one hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of (True, reason) if classified as vitamin D, else (False, reason).
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (250-600 Da is typical for vitamin D)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.1f} Da is out of vitamin D range (250-600 Da)"
    
    # Look for at least one hydroxyl (-OH) group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl (-OH) group found"
    
    # Get ring information and all rings as lists of atom indices
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Build a graph of rings: node for each ring; two rings are connected if they share >=2 atoms.
    n = len(all_rings)
    adjacency = {i: set() for i in range(n)}
    ring_sets = [set(r) for r in all_rings]
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Find connected components (clusters) of rings using DFS
    visited = set()
    clusters = []
    for i in range(n):
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
    
    # Get the largest ring cluster
    largest_cluster = max(clusters, key=lambda comp: len(comp))
    cluster_ring_count = len(largest_cluster)
    
    # Accept clusters with 2 to 5 rings (allowing for secosteroid cores that lost one ring)
    if cluster_ring_count < 2 or cluster_ring_count > 5:
        return False, f"Fused ring cluster has {cluster_ring_count} rings; expected between 2 and 5 for a vitamin D secosteroid core"
    
    # Gather all atom indices in the largest cluster
    cluster_atoms = set()
    for idx in largest_cluster:
        cluster_atoms.update(ring_sets[idx])
    
    # Check aromaticity of the fused cluster. Vitamin D cores are mostly non‐aromatic.
    arom_count = 0
    for atom_idx in cluster_atoms:
        if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
            arom_count += 1
    aromatic_ratio = arom_count / len(cluster_atoms)
    if aromatic_ratio > 0.5:
        return False, f"Fused ring cluster is predominantly aromatic (ratio {aromatic_ratio:.2f}); not consistent with vitamin D secosteroid core"
    
    # Now search for an extended aliphatic side chain attached to the fused cluster.
    # We require a branch (found on an atom in the cluster) that is at least 3 carbon atoms long
    # and has at least one hydroxyl group (that is not part of the fused cluster).
    def branch_has_oh_and_min_c(start_idx, min_c=3, max_depth=8):
        """Explore branch starting from start_idx (which is not in the cluster)
           and coming off a cluster atom. Return True if a chain with at least min_c carbon atoms 
           and containing an -OH is found.
        """
        visited_atoms = set()
        stack = [(start_idx, 0, 0)]  # (atom_index, current depth, carbon_count)
        while stack:
            atom_idx, depth, c_count = stack.pop()
            if depth > max_depth:
                continue
            if atom_idx in visited_atoms:
                continue
            visited_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            # Increment carbon count if carbon
            if atom.GetAtomicNum() == 6:
                c_count += 1
            # Check if this atom is a hydroxyl oxygen (with H)
            if atom.GetAtomicNum() == 8:
                # A simple check for -OH: oxygen attached to at least one hydrogen.
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 1:
                        if c_count >= min_c:
                            return True
            # Continue traversal: avoid jumping into any atom already in the fused cluster.
            for nb in atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                if nb_idx in cluster_atoms:
                    continue
                stack.append((nb_idx, depth+1, c_count))
        return False

    branch_found = False
    # For each atom in the fused cluster, check neighbors not in the cluster.
    for atom_idx in cluster_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in cluster_atoms:
                continue
            # Only search branch if the neighbor is aliphatic (carbon or oxygen is acceptable)
            if nb.GetAtomicNum() not in (6, 8):
                continue
            if branch_has_oh_and_min_c(nb_idx):
                branch_found = True
                break
        if branch_found:
            break
    if not branch_found:
        return False, "No extended aliphatic side chain (with at least 3 carbons and a hydroxyl) attached to the fused ring cluster"
    
    # If all tests pass, classify as a vitamin D compound.
    return True, "Matches vitamin D criteria: fused non‐aromatic secosteroid core with hydroxyl group(s) and appropriate side chain"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Testing with one example: 24-fluoro-25-hydroxyvitamin D3
    test_smiles = "FC(CC[C@H]([C@@H]1[C@@]2([C@@H](CC1)/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)C)C)C(O)(C)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)
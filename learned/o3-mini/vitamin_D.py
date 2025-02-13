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

Heuristic criteria used here:
 1. Valid molecule from SMILES.
 2. Has at least one hydroxyl (-OH) group.
 3. Contains a fused ring system resembling a secosteroid core.
    – Instead of forcing exactly 3 rings we look for a “cluster” of 3–5 rings that share bonds.
    – We further require that the fused system is not heavily aromatic (vitamin D cores are not fully aromatic).
 4. Molecular weight is in the expected range (250–600 Da).
 
Note:
This heuristic is an approximation; many vitamin D compounds may have slight variations and extra ring(s)
(e.g. from epoxides or lactol rings) and many non‐vitamin D compounds (such as fully aromatic 3‐ring systems)
may falsely match a simpler rule. We attempt to reduce false positives by combining ring‐fusion analysis
with an “aromaticity” filter.
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
      - Contains a fused ring cluster (a group of rings sharing edges) where the number of rings is 
        typically 3–5 and the rings are not mostly aromatic. (Steroid cores are saturated or partially unsaturated.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of (True, reason) if the molecule is classified as a vitamin D compound,
                     otherwise (False, reason) explaining why it was rejected.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.1f} Da is out of expected vitamin D range (250-600 Da)"
    
    # Check for at least one hydroxyl (–OH) group using a SMARTS pattern
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl (-OH) group found"
    
    # Identify all rings in the molecule
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Build a graph of rings: each ring is a node; add an edge if two rings share at least two atoms
    n = len(all_rings)
    adjacency = {i: set() for i in range(n)}
    ring_sets = [set(r) for r in all_rings]
    for i in range(n):
        for j in range(i + 1, n):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                adjacency[i].add(j)
                adjacency[j].add(i)
                
    # Find connected components (clusters) of rings by DFS
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
            
    # Get the largest ring cluster (by number of rings)
    largest_cluster = max(clusters, key=lambda comp: len(comp))
    cluster_ring_count = len(largest_cluster)
    
    # We expect a secosteroid core to have a fused ring cluster of around 3-5 rings.
    if cluster_ring_count < 3 or cluster_ring_count > 5:
        return False, f"Fused ring cluster has {cluster_ring_count} rings; expected between 3 and 5 for a secosteroid core"
    
    # Collect all atom indices in the largest cluster (union over each ring in the cluster)
    cluster_atoms = set()
    for idx in largest_cluster:
        cluster_atoms.update(ring_sets[idx])
    
    # Evaluate the “aromaticity” of the cluster. In steroids the rings are mostly aliphatic.
    arom_count = 0
    for atom_idx in cluster_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic():
            arom_count += 1
    aromatic_ratio = arom_count / len(cluster_atoms)
    # If most of the atoms in the fused cluster are aromatic, that is not steroid-like.
    if aromatic_ratio > 0.5:
        return False, f"Fused ring cluster is predominantly aromatic (aromatic ratio {aromatic_ratio:.2f}); not consistent with vitamin D secosteroid core"
    
    # Passed all checks: the molecule is a plausible vitamin D compound.
    return True, "Matches vitamin D criteria: fused non‐aromatic secosteroid core with hydroxyl group(s) and proper molecular weight"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Testing with one example: 24-fluoro-25-hydroxyvitamin D3
    test_smiles = "FC(CC[C@H]([C@@H]1[C@@]2([C@@H](CC1)/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)C)C)C(O)(C)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)
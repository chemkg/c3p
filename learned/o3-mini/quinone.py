"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: Quinones, defined as compounds having a fully conjugated cyclic dione structure.
For example, benzoquinones and polycyclic/heterocyclic analogues.
This heuristic attempts to merge individual SSSR rings into fused ring systems,
then inspects each system for:
  - A minimum total size (≥ 5 atoms)
  - A high fraction (≥90%) of ring atoms that are aromatic or sp2‐hybridized
  - Exactly two carbonyl groups (C(=O)) embedded in the system. In this context,
    a carbonyl group is counted only if the carbon is part of the ring system and
    has a double bond to an oxygen that is outside.
  - Each carbonyl carbon is required to have at least 2 neighboring atoms in the system
    that are aromatic or sp2‐hybridized.
  - The two carbonyl atoms are not “too close” in the fused system; the shortest bond path 
    (restricted to atoms of the fused system) is at least 3.
Note: Given the wide range of quinone structures, this remains only a heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdchem
from collections import deque

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    The heuristic builds fused ring systems from the SSSR rings and then for each system inspects
    the level of conjugation and verifies the presence of exactly 2 suitably separated carbonyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule classified as a quinone, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass  # proceed even if sanitization fails

    # Get rings from RDKit (each as a tuple of atom indices)
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())
    if not rings:
        return False, "No ring systems found"

    # Merge rings that have any overlapping atoms to get fused ring systems.
    # We will iteratively merge rings that share atoms.
    def merge_rings(rings_list):
        merged = [set(r) for r in rings_list]
        merged_changed = True
        while merged_changed:
            merged_changed = False
            new_merged = []
            while merged:
                current = merged.pop(0)
                merge_occurred = False
                # Look for any other set that overlaps current
                i = 0
                while i < len(merged):
                    other = merged[i]
                    if current & other:
                        # Merge and remove the other set
                        current |= other
                        merged.pop(i)
                        merge_occurred = True
                        merged_changed = True
                    else:
                        i += 1
                new_merged.append(current)
            merged = new_merged
        return merged

    fused_systems = merge_rings(rings)

    # Helper function: check if an atom (assumed to be carbon) is a carbonyl in a given ring system.
    # It must be in the ring set and have a double bond to oxygen that is NOT in the ring.
    def is_ring_carbonyl(atom, ring_set):
        if atom.GetAtomicNum() != 6:
            return False
        has_dblO = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    # Only count if the oxygen is not part of the ring.
                    if nbr.GetIdx() not in ring_set:
                        has_dblO = True
        return has_dblO

    # For each fused ring system, check our criteria.
    for ring_set in fused_systems:
        if len(ring_set) < 5:
            continue

        # Get the atoms in the fused system.
        system_atoms = [mol.GetAtomWithIdx(idx) for idx in ring_set]
        # Count atoms that are aromatic or sp2 hybridized.
        n_conj = 0
        for atom in system_atoms:
            if atom.GetIsAromatic() or atom.GetHybridization() == rdchem.HybridizationType.SP2:
                n_conj += 1
        if n_conj < 0.9 * len(system_atoms):
            continue

        # Identify carbonyl atoms within this fused system.
        carbonyl_atoms = []
        for atom in system_atoms:
            if is_ring_carbonyl(atom, ring_set):
                carbonyl_atoms.append(atom)

        if len(carbonyl_atoms) != 2:
            continue

        # For each carbonyl atom, check that at least 2 of its neighbors inside the fused system
        # are conjugated (i.e. aromatic or sp2).
        valid = True
        for atom in carbonyl_atoms:
            in_system_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_set]
            n_conj_neighbors = sum(1 for nbr in in_system_neighbors if nbr.GetIsAromatic() or nbr.GetHybridization() == rdchem.HybridizationType.SP2)
            if n_conj_neighbors < 2:
                valid = False
                break
        if not valid:
            continue

        # Compute bond distance between the two carbonyl atoms, restricted to the fused system.
        # Build a graph (adjacency dict) for atoms in the ring system.
        graph = {idx: set() for idx in ring_set}
        for idx in ring_set:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring_set:
                    graph[idx].add(nbr_idx)

        start = carbonyl_atoms[0].GetIdx()
        end = carbonyl_atoms[1].GetIdx()
        # Use breadth-first search (BFS) to compute shortest path length.
        visited = {start}
        queue = deque([(start, 0)])
        found_path_length = None
        while queue:
            current, dist = queue.popleft()
            if current == end:
                found_path_length = dist
                break
            for nbr in graph[current]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, dist + 1))
        if found_path_length is None or found_path_length < 3:
            # The carbonyls are too close or no path found.
            continue

        # If we reach here, the fused ring system qualifies as a quinone.
        return True, "Found a fused ring system with two carbonyl groups embedded in a fully conjugated cyclic system"

    return False, "No fully conjugated cyclic dione (quinone) structure was found"


# For testing purposes (optional)
if __name__ == "__main__":
    test_smiles = [
        "COC1=C(OC)C(=O)C(C)=CC1=O",   # ubiquinone-0 (expected True)
        "O=C1C=C(OC)[C@](O)(C)C([C@@]1(OC(=O)C(CC)C)C)=O",  # Phomaligol A (should be classified as quinone)
        "C1(C=C(C(C(=C1)Cl)=O)O)=O",    # 2-chloro-6-hydroxy-1,4-benzoquinone (expected True)
        "O=C1N(NC(=O)C=C1C)C2=CC=CC=C2", # example false positive in previous run
    ]
    for s in test_smiles:
        result, reason = is_quinone(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*40}")
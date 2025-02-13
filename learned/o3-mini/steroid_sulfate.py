"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
The molecule must contain a steroid nucleus (a fused tetracyclic system where the fused core is composed mostly of carbon) 
with a sulfate group attached via an oxygen to one of the steroid carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    The following conditions must be met:
      1. The molecule must contain a steroid nucleus.
         Here, we approximate this feature by finding a fused ring cluster that contains at least four rings.
         We then require that the majority of the atoms in the combined fused rings are carbons.
      2. The molecule must contain a sulfate ester group (OS(=O)(=O) or its ionic form) that is attached via its oxygen
         to one of the atoms in the steroid nucleus (indicating derivatization of a hydroxy group).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a steroid sulfate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- STEP 1: Identify a steroid nucleus by analyzing fused rings ---
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices in each ring
    if not atom_rings:
        return False, "No rings found in the molecule; cannot be a steroid."
    
    # Build a graph where each node is a ring (index in atom_rings) and edges connect rings that are fused.
    # Two rings are considered fused if they share at least two atoms (i.e. at least one bond).
    ring_graph = {i: set() for i in range(len(atom_rings))}
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            # Check intersection of ring i and ring j
            if len(set(atom_rings[i]).intersection(atom_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components (clusters) of fused rings.
    visited = set()
    fused_clusters = []
    for i in range(len(atom_rings)):
        if i not in visited:
            stack = [i]
            cluster = set()
            while stack:
                cur = stack.pop()
                if cur not in visited:
                    visited.add(cur)
                    cluster.add(cur)
                    stack.extend(ring_graph[cur] - visited)
            fused_clusters.append(cluster)
    
    # Look for any cluster with at least 4 fused rings and validate its atom composition.
    steroid_core_atoms = None
    for cluster in fused_clusters:
        if len(cluster) >= 4:
            # Get union of all atoms in rings of this cluster.
            atoms_in_cluster = set()
            for ring_idx in cluster:
                atoms_in_cluster.update(atom_rings[ring_idx])
            # Count carbon atoms in this cluster
            num_carbons = sum(1 for idx in atoms_in_cluster if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            ratio_carbons = num_carbons / len(atoms_in_cluster)
            # For a steroid, we expect most atoms in the fused rings to be carbon (say at least 70%).
            if ratio_carbons >= 0.7:
                steroid_core_atoms = atoms_in_cluster
                break
    if steroid_core_atoms is None:
        return False, "No steroid-like fused ring system (four fused rings with mostly carbons) found."
        
    # --- STEP 2: Identify a sulfate ester group attached to the steroid nucleus ---
    # Define a SMARTS pattern for sulfate ester group.
    # This pattern will match an oxygen single-bonded to sulfur with two double bonds to oxygen.
    # It should match both neutral and ionic forms.
    sulfate_pattern = Chem.MolFromSmarts("O[S;$([OX2])](=O)(=O)")  # The [OX2] ensures the oxygen has two connections.
    if sulfate_pattern is None:
        return False, "Error in sulfate SMARTS pattern."
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate ester group found."
    
    # For each sulfate match, check that its oxygen (first atom in match) is attached to a carbon
    # that is part of the steroid nucleus.
    for match in sulfate_matches:
        sulfate_oxygen = mol.GetAtomWithIdx(match[0])
        # Examine neighbors (skip the S atom that is part of the sulfate)
        for nbr in sulfate_oxygen.GetNeighbors():
            # If neighbor is sulfur (atomic num 16) then skip it
            if nbr.GetAtomicNum() == 16:
                continue
            # If neighbor is a carbon and is in the steroid nucleus, assume it was originally a hydroxy on the steroid.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in steroid_core_atoms:
                return True, "Molecule contains a steroid nucleus with a sulfate ester derived from a steroid hydroxy group."
    
    return False, "Sulfate ester group not attached to the steroid nucleus as required."

# For testing purposes (you can run this to perform a simple test):
if __name__ == "__main__":
    test_smiles = "O(S(O)(=O)=O)[C@@H]1CC=2[C@]([C@]3(CC[C@]4([C@]([C@@]3(CC2)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])(C)CC1"
    result, reason = is_steroid_sulfate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
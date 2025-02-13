"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: Glycosaminoglycan
Definition: Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.
Heuristic:
    1. Identify candidate sugar rings defined as 5- or 6-membered rings that contain exactly one ring oxygen.
    2. For each candidate ring, check for an external amino substituent (an atom with atomic number 7 not in the ring).
    3. Build a connectivity graph between candidate sugar rings: two rings are “linked” if any atom in one is directly bonded
       to any atom in the other.
    4. Only if there is a sufficiently large (>= 3 rings) connected polysaccharide chain AND
       at least 30% of the candidate rings show an external amino substituent do we classify the molecule as a glycosaminoglycan.
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.

    Glycosaminoglycans are defined as polysaccharides (a chain of three or more carbohydrate units)
    that contain a substantial proportion of aminomonosaccharide residues.
    
    The heuristic is as follows:
     - Identify candidate sugar rings: rings of size 5 or 6 that have exactly one oxygen.
     - Flag a candidate sugar ring as “aminosugar” if at least one of its ring atoms has a substituent (not in the ring)
       that is a nitrogen atom.
     - Build a connectivity graph between rings: two rings are connected if any atom in one ring is directly bonded
       to any atom in the other.
     - Require that the molecule has at least 2 candidate sugar rings overall AND that at least one connected cluster
       (as measured by the size of the connected component in the graph) is of size 3 or more (i.e. a polysaccharide chain).
     - Finally, require that at least 30% of all candidate sugar rings carry an amino substituent.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a glycosaminoglycan, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []       # List of candidate rings (each is a set of atom indices)
    aminosugar_flags = []      # List of booleans: True if the candidate ring has an external amino group

    # Loop over all rings in the molecule.
    for ring in atom_rings:
        # Only consider rings of size 5 or 6.
        if len(ring) not in (5,6):
            continue

        # Count the number of oxygen atoms in this ring.
        oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
        # Many carbohydrate rings (furanoses or pyranoses) have exactly one ring oxygen.
        if oxygen_count != 1:
            continue

        # Add this ring (as a set of atom indices) to our candidate list.
        ring_set = set(ring)
        candidate_rings.append(ring_set)
        
        # Check for an external amino substituent on any atom of the ring.
        found_amino = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Inspect neighbors outside the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetAtomicNum() == 7:
                    found_amino = True
                    break
            if found_amino:
                break
        aminosugar_flags.append(found_amino)
    
    total_rings = len(candidate_rings)
    if total_rings < 2:
        return False, f"Only {total_rings} candidate sugar ring(s) found; need a polysaccharide chain."
    
    # Build a connectivity graph between candidate rings.
    # Two candidate rings are connected if any atom in one ring is directly bonded to any atom in the other.
    n = total_rings
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # Check if any atom in candidate_rings[i] is directly bonded to any atom in candidate_rings[j]
            connected = False
            for idx in candidate_rings[i]:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in candidate_rings[j]:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                graph[i].add(j)
                graph[j].add(i)
    
    # Compute connected components in the graph.
    visited = set()
    components = []
    for i in range(n):
        if i in visited:
            continue
        stack = [i]
        comp = set()
        while stack:
            node = stack.pop()
            if node in comp:
                continue
            comp.add(node)
            for neighbor in graph[node]:
                if neighbor not in comp:
                    stack.append(neighbor)
        visited |= comp
        components.append(comp)
    
    largest_comp_size = max(len(comp) for comp in components) if components else 0
    if largest_comp_size < 3:
        return False, (f"Largest connected candidate sugar ring component has {largest_comp_size} ring(s); "
                       "need at least 3 connected rings to indicate a polysaccharide chain.")
    
    # Now compute the fraction of candidate rings that are aminosugars.
    count_amino = sum(1 for flag in aminosugar_flags if flag)
    fraction_amino = count_amino / total_rings
    if fraction_amino < 0.3:
        return False, (f"Only {count_amino} out of {total_rings} candidate sugar rings "
                       f"({fraction_amino*100:.1f}%) have an amino substituent; "
                       "not enough aminomonosaccharide residues.")
    
    return True, (f"Detected {total_rings} candidate sugar rings, with {count_amino} "
                  f"aminosugar(s) ({fraction_amino*100:.1f}%), and a largest connected chain of {largest_comp_size} rings, "
                  "consistent with a glycosaminoglycan.")

# Example usage:
if __name__ == "__main__":
    # A test SMILES of a glycosaminoglycan fragment (this is only an example and may not cover full complexity).
    test_smiles = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(CO)C(O)C(O)C2O)[C@@H](O)[C@@H]1O"
    result, reason = is_glycosaminoglycan(test_smiles)
    print(result, reason)
"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: Glycosaminoglycan
Definition: Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.
Heuristic improvements:
    1. Identify candidate sugar rings as 5‐ or 6‐membered rings that contain exactly one ring oxygen.
    2. For a ring to be considered “sugar‐like” it must have at least one external substituent 
       that is either a hydroxyl (-OH) (oxygen with at least one hydrogen) or an amino group (-NH2) 
       (nitrogen with at least one attached hydrogen).
    3. Flag candidate rings as aminomonosaccharides if (among the substituents) at least one is an –NH2.
    4. Build a connectivity graph among candidate rings by requiring that two rings are “linked” 
       if one of the ring atoms is connected to an external oxygen (i.e. a potential glycosidic linker) 
       that in turn is connected to an atom in the other ring OR the rings are directly bonded.
    5. Require that the largest connected component (polysaccharide chain) has at least 4 rings, and 
       that ≥30% of candidate rings exhibit an amino substituent.
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
    Heuristic:
      - First, we identify candidate sugar rings as 5- or 6-membered rings that have exactly one oxygen atom.
      - Additionally, at least one ring atom must have an external substituent that is either a hydroxyl (-OH)
        (oxygen with at least one hydrogen) or an amino group (-NH2; nitrogen with at least one hydrogen).
      - The ring is flagged as an aminomonosaccharide if one or more of these substituents is nitrogen.
      - Two candidate rings are considered connected if either they are directly bonded OR 
        they are linked by an oxygen (i.e. one of the ring atoms is bonded to an external oxygen which then is 
        bonded to another candidate ring). This mimics the glycosidic bond.
      - Finally, we require that:
            (a) the molecule has at least 2 candidate sugar rings,
            (b) the largest connected cluster of candidate rings has at least 4 rings,
            (c) at least 30% of candidate rings are aminomonosaccharides.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a glycosaminoglycan, False otherwise.
      str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve all ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []      # each candidate ring: a set of atom indices of the ring
    amino_flags = []          # parallel list: True if ring has an external amino substituent
    
    # Loop over all rings in the molecule.
    for ring in atom_rings:
        # Consider only rings of size 5 or 6.
        if len(ring) not in (5,6):
            continue
        
        # Check that the ring has exactly one oxygen atom.
        oxygen_in_ring = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_in_ring += 1
        if oxygen_in_ring != 1:
            continue
        
        # Look for external substituents on ring atoms.
        # Expect at least one substituent to be an -OH (oxygen with H) or -NH2 (nitrogen with H).
        has_valid_substituent = False
        has_amino = False
        ring_set = set(ring)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors not in the ring (external substituents)
                if nbr.GetIdx() in ring_set:
                    continue
                # For oxygen substituents: check if it likely represents a hydroxyl: has one or more H.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    has_valid_substituent = True
                # For nitrogen: check for at least one hydrogen (to distinguish from e.g. an amide N)
                if nbr.GetAtomicNum() == 7 and nbr.GetTotalNumHs() > 0:
                    has_valid_substituent = True
                    has_amino = True
            # If we already found a valid substituent, then no need to check further for this ring.
            if has_valid_substituent:
                break
        
        # If the ring does not show any hydroxyl or amino substituent, discard it.
        if not has_valid_substituent:
            continue
        
        candidate_rings.append(ring_set)
        amino_flags.append(has_amino)
    
    total_rings = len(candidate_rings)
    if total_rings < 2:
        return False, f"Only {total_rings} candidate sugar ring(s) found; need a longer polysaccharide chain."
    
    # Build connectivity graph between candidate rings.
    # We add an edge if the candidate rings are directly bonded or if they are bridged by an oxygen.
    n = total_rings
    graph = {i: set() for i in range(n)}
    
    # First, try direct connection (if any atom in ring i is bonded directly to an atom in ring j).
    for i in range(n):
        for j in range(i+1, n):
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
    
    # Now, also add connectivity via a bridging oxygen.
    # For each candidate ring, check each ring atom's external neighbors that are oxygen.
    # Then if that oxygen has another neighbor belonging to another candidate ring, mark connection.
    for i in range(n):
        for idx in candidate_rings[i]:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in candidate_rings[i]:
                    continue  # skip atoms within same ring
                # Consider bridging only if the neighbor is oxygen.
                if nbr.GetAtomicNum() == 8:
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == atom.GetIdx():
                            continue
                        # Check if this neighbor belongs to any candidate ring j (different from i).
                        for j in range(n):
                            if j == i:
                                continue
                            if nbr2.GetIdx() in candidate_rings[j]:
                                # Optionally, one could require that the bridging oxygen is not a hydroxyl on a substituent
                                # but rather part of an ether linkage. Here we assume any oxygen bridging two rings qualifies.
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
    
    largest_component = max(components, key=len) if components else set()
    largest_size = len(largest_component)
    if largest_size < 4:
        return False, (f"Largest connected candidate sugar ring component has {largest_size} rings; "
                       "need at least 4 connected rings to indicate a glycosaminoglycan.")
    
    # Calculate the fraction of candidate rings that are aminomonosaccharides.
    count_amino = sum(1 for flag in amino_flags if flag)
    fraction_amino = count_amino / total_rings
    if fraction_amino < 0.3:
        return False, (f"Only {count_amino} out of {total_rings} candidate sugar rings "
                       f"({fraction_amino*100:.1f}%) have an amino substituent; "
                       "not enough aminomonosaccharide residues.")
    
    return True, (f"Detected {total_rings} candidate sugar rings, with {count_amino} "
                  f"aminosugar(s) ({fraction_amino*100:.1f}%), and a largest connected chain of {largest_size} rings, "
                  "consistent with a glycosaminoglycan.")

# Example usage: (this block can be used for testing)
if __name__ == "__main__":
    # Test with a putative glycosaminoglycan fragment (example SMILES, may not be comprehensive)
    test_smiles = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(CO)C(O)C(O)C2O)[C@@H](O)[C@@H]1O"
    result, reason = is_glycosaminoglycan(test_smiles)
    print(result, reason)
"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: Glycosaminoglycan
Definition: Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.
Heuristic improvements in this version:
  1. Identify candidate sugar rings as 5‐ or 6‐membered rings with exactly one ring oxygen.
  2. Each candidate ring must have at least one external substituent (attached via a single bond)
     that is either a hydroxyl group (oxygen with at least one hydrogen) or an amino group (nitrogen with at least one hydrogen).
  3. Flag candidate rings as aminomonosaccharide if any valid substituent is nitrogen.
  4. When building a connectivity graph among candidate rings, add an edge if:
       (a) any atom in one ring is directly bonded to an atom in another ring, OR
       (b) there is an external oxygen (with degree exactly 2) attached via a single bond from a ring atom that bridges to an atom in a different candidate ring.
  5. Require that the molecule has at least 2 candidate rings,
     the largest connected cluster has at least 4 rings, and
     ≥30% of candidate rings show an amino substituent.
  
If determining these features is too challenging, the function may return (None, None).
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
    Returns:
      bool: True if classified as a glycosaminoglycan, False otherwise.
      str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []  # list of sets of atom indices representing candidate sugar rings
    amino_flags = []      # parallel list: True if candidate ring has an -NH2 substituent
    
    for ring in atom_rings:
        # Consider only rings of size 5 or 6.
        if len(ring) not in (5,6):
            continue
        ring_set = set(ring)
        # Count oxygen atoms in the ring.
        ring_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if ring_oxygens != 1:
            continue
        
        # Look for external substituents (attached by a single bond).
        valid_substituent_found = False
        amino_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors not in the ring and attached via a single bond
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetIdx() in ring_set or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # Check for hydroxyl: oxygen with at least one hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    valid_substituent_found = True
                # Check for amino: nitrogen with at least one hydrogen.
                if nbr.GetAtomicNum() == 7 and nbr.GetTotalNumHs() > 0:
                    valid_substituent_found = True
                    amino_found = True
            # If we already found a valid substituent on one atom, we can stop checking further in this ring.
            if valid_substituent_found:
                break

        if not valid_substituent_found:
            continue
        
        candidate_rings.append(ring_set)
        amino_flags.append(amino_found)
    
    total_candidates = len(candidate_rings)
    if total_candidates < 2:
        return False, f"Only {total_candidates} candidate sugar ring(s) found; need a longer polysaccharide chain."
    
    # Build graph connectivity among candidate rings.
    # Nodes are indexes into candidate_rings.
    n = total_candidates
    graph = {i: set() for i in range(n)}
    
    # Add edge if rings share a bond (direct connection).
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
    
    # Now add connectivity via bridging oxygen.
    # For each candidate ring, for each external neighbor that is oxygen and attached via a single bond,
    # if that oxygen atom has degree exactly 2 and its other neighbor belongs to a different candidate ring, add an edge.
    for i in range(n):
        for idx in candidate_rings[i]:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in candidate_rings[i]:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                if nbr.GetAtomicNum() != 8:
                    continue
                if nbr.GetDegree() != 2:
                    continue
                # Get the other neighbor of the bridging oxygen.
                nbr_atoms = [n_atom for n_atom in nbr.GetNeighbors() if n_atom.GetIdx() != atom.GetIdx()]
                if not nbr_atoms:
                    continue
                other = nbr_atoms[0]
                # Check for every other candidate ring j if the other neighbor belongs to ring j.
                for j in range(n):
                    if j == i:
                        continue
                    if other.GetIdx() in candidate_rings[j]:
                        graph[i].add(j)
                        graph[j].add(i)
    
    # Compute connected components of the candidate ring graph.
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
    
    if not components:
        return False, "No connected candidate ring clusters found."
    
    largest_component = max(components, key=len)
    largest_size = len(largest_component)
    if largest_size < 4:
        return False, f"Largest connected candidate sugar ring component has {largest_size} rings; need at least 4 connected rings to indicate a glycosaminoglycan."
    
    count_amino = sum(1 for flag in amino_flags if flag)
    fraction_amino = count_amino / total_candidates
    if fraction_amino < 0.3:
        return False, (f"Only {count_amino} out of {total_candidates} candidate sugar rings "
                       f"({fraction_amino * 100:.1f}%) have an amino substituent; not enough aminomonosaccharide residues.")
    
    return True, (f"Detected {total_candidates} candidate sugar rings, with {count_amino} aminosugar(s) "
                  f"({fraction_amino * 100:.1f}%), and a largest connected chain of {largest_size} rings, "
                  "consistent with a glycosaminoglycan.")

# For testing purposes:
if __name__ == "__main__":
    # Example SMILES (this test fragment is only for demonstration and may need adjustment)
    test_smiles = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2OC(CO)C(O)C(O)C2O)[C@@H](O)[C@@H]1O"
    result, reason = is_glycosaminoglycan(test_smiles)
    print(result, reason)
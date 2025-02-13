"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition: An oligosaccharide comprising four monomeric monosaccharide units.
This revised classifier:
  - Parses the SMILES and requires one connected component.
  - Detects candidate sugar rings as rings of size 5 (furanose) or 6 (pyranose) containing exactly one oxygen.
  - Requires that each candidate sugar ring has at least one exocyclic oxygen substituent carrying a hydrogen (e.g. –OH)
    and at least one chiral atom.
  - In addition, when checking connectivity, the candidate ring’s set is “grown” by including neighboring oxygen atoms
    that commonly act as glycosidic linkers.
If there are exactly 4 candidate rings and they are connected via such linkers, the molecule qualifies as a tetrasaccharide.
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determine if a molecule (given as a SMILES string) is a tetrasaccharide.
    The procedure:
      - Ensures the SMILES is valid and corresponds to a single (connected) molecule.
      - From the ring information, it gathers candidate sugar rings defined as:
          * A ring of 5 or 6 atoms,
          * Containing exactly one oxygen atom (and the remaining being carbons),
          * Having at least one exocyclic –OH substituent (an oxygen neighbor having at least one hydrogen),
          * And having at least one explicitly chiral center.
      - It then “grows” each candidate unit’s atom set by including neighbouring oxygen atoms (the typical glycosidic linker atoms).
      - Finally, it requires exactly 4 candidate sugar rings whose grown sets are mutually connected (i.e. intersect at least one other).
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): (True, reason) if the molecule qualifies as a tetrasaccharide;
                   (False, reason) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a single connected component
    mol_smiles = Chem.MolToSmiles(mol)
    if '.' in mol_smiles:
        return False, "Molecule contains multiple disconnected fragments"
    
    ring_info = mol.GetRingInfo().AtomRings()
    
    def is_candidate_sugar_ring(atom_indices):
        """
        Evaluate whether a given ring (identified by its atom indices) meets the criteria for a
        candidate sugar ring.
          - Ring must be of size 5 (furanose) or 6 (pyranose).
          - The ring must have exactly one oxygen (rest assumed to be carbons).
          - At least one exocyclic -OH substituent: one neighbor (not in the ring) that is an oxygen
            with at least one hydrogen.
          - At least one chiral center should be specified in the ring.
        """
        n = len(atom_indices)
        if n not in (5, 6):
            return False
        oxy_count = 0
        carbon_count = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxy_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_count += 1
        # Check typical sugar ring composition: one oxygen and the rest carbons.
        if oxy_count != 1 or (n - oxy_count) != carbon_count:
            return False
        
        # Check for exocyclic hydroxyl group(s)
        exo_oh = 0
        for idx in atom_indices:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in atom_indices:
                    continue
                # We expect an -OH substituent: neighbor oxygen with at least one hydrogen
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    exo_oh += 1
        if exo_oh < 1:
            return False
        
        # Relaxed chiral requirement: at least one chiral center in the ring.
        chiral = 0
        for idx in atom_indices:
            a = mol.GetAtomWithIdx(idx)
            if a.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                chiral += 1
        if chiral < 1:
            return False
        
        return True

    # Gather candidate sugar rings (as sets of atom indices)
    candidate_rings = []
    for ring in ring_info:
        if is_candidate_sugar_ring(ring):
            candidate_rings.append(set(ring))
            
    if len(candidate_rings) != 4:
        return False, f"Found {len(candidate_rings)} candidate carbohydrate rings, need exactly 4"

    # Build a connectivity graph among candidate rings.
    # Instead of checking only for direct bonds between atoms in rings, we "grow" the candidate set by including
    # neighboring oxygen atoms (which may serve as glycosidic linkers).
    grown_sets = []
    for ring_atoms in candidate_rings:
        grown = set(ring_atoms)  # start with atoms in the ring
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # If neighbor is oxygen, add it to the grown set.
                if nbr.GetAtomicNum() == 8:
                    grown.add(nbr.GetIdx())
        grown_sets.append(grown)
    
    n = len(grown_sets)
    connectivity_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # If the grown sets of candidate rings share any atom, treat them as connected.
            if grown_sets[i].intersection(grown_sets[j]):
                connectivity_graph[i].add(j)
                connectivity_graph[j].add(i)
    
    # Do a DFS to check if all four candidate rings are connected.
    visited = set()
    def dfs(node):
        stack = [node]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            stack.extend(list(connectivity_graph[cur] - visited))
    
    dfs(0)
    if len(visited) != n:
        return False, "Candidate carbohydrate rings are not connected via glycosidic bonds"
    
    return True, "Molecule contains exactly 4 connected monosaccharide (sugar) units"

# For testing purposes (script usage)
if __name__ == '__main__':
    # A known tetrasaccharide example as provided.
    test_smiles = "OC[C@H]1O[C@@H](OC[C@H]2O[C@@H](OC[C@H]3O[C@@H](OC[C@H]4OC(O)[C@H](O)[C@@H](O)[C@H]4O)[C@H](O)[C@@H](O)[C@H]3O)[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
    result, reason = is_tetrasaccharide(test_smiles)
    print(result, reason)
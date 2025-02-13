"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition: An oligosaccharide comprising four monomeric monosaccharide units.
This revised classifier:
  - Requires a valid single-component molecule.
  - Searches the molecule’s rings for candidate sugar rings.
    For our purposes a candidate sugar ring is defined as:
      * A 5- or 6-membered ring,
      * Containing exactly one oxygen atom (the ring heteroatom) and otherwise carbons,
      * Having at least one exocyclic hydroxyl group (an oxygen neighbor outside the ring that bears at least one hydrogen).
      * (We relax requirements on chiral specification.)
  - Each candidate ring is “grown” into a unit by including all neighboring oxygen atoms (the typical glycosidic linkers).
  - Finally, if exactly 4 candidate units are found and the “grown” sets are mutually connected (via intersection of their atoms), 
    the molecule qualifies as a tetrasaccharide.
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determine if a molecule (given as a SMILES string) is a tetrasaccharide.
    The procedure:
      - Parse SMILES and ensure the molecule is a single (connected) component.
      - From the ring information identify candidate sugar rings defined as:
           • 5- or 6-membered ring
           • Exactly one ring oxygen (atomic number 8) and the rest carbons
           • At least one exocyclic –OH group (an external oxygen that has at least one hydrogen)
         (Chirality is not required explicitly.)
      - “Grow” each candidate ring by adding any immediate neighboring oxygens.
      - Build a connectivity graph among the four candidate units; if they are connected via shared (or neighboring) oxygen atoms,
        the molecule qualifies.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): (True, reason) if the molecule qualifies as a tetrasaccharide; otherwise (False, reason).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is a single connected component.
    mol_smiles = Chem.MolToSmiles(mol)
    if '.' in mol_smiles:
        return False, "Molecule contains multiple disconnected fragments"
    
    ring_info = mol.GetRingInfo().AtomRings()

    def is_candidate_sugar_ring(atom_indices):
        """
        Evaluate whether a given ring (identified by its atom indices) meets the criteria for a candidate sugar ring.
          - Ring size must be 5 (furanose) or 6 (pyranose).
          - Among the ring atoms, exactly one must be oxygen (atomic number 8) and the remainder carbon (atomic number 6).
          - At least one exocyclic –OH group must be attached to an atom in the ring;
            that is, a neighbor (not in the ring) that is oxygen and carries at least one hydrogen.
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
            else:
                # If any atom in the ring is neither carbon nor oxygen,
                # then it does not meet our idealized sugar ring criteria.
                return False
        if oxy_count != 1 or (n - oxy_count) != carbon_count:
            return False
        
        # Check for exocyclic hydroxyl group(s):
        exo_oh_found = False
        for idx in atom_indices:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in atom_indices:
                    continue
                # We require an -OH group, i.e. a neighbor that is oxygen and has at least one hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    exo_oh_found = True
                    break
            if exo_oh_found:
                break
        if not exo_oh_found:
            return False

        return True

    # Gather candidate sugar rings (each as a set of atom indices).
    candidate_rings = []
    for ring in ring_info:
        if is_candidate_sugar_ring(ring):
            candidate_rings.append(set(ring))
            
    if len(candidate_rings) != 4:
        return False, f"Found {len(candidate_rings)} candidate carbohydrate rings, need exactly 4"
    
    # For each candidate ring, “grow” its set by adding any neighboring oxygen atoms.
    grown_sets = []
    for ring_atoms in candidate_rings:
        grown = set(ring_atoms)  # start with the ring atoms
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # If neighbor is oxygen (and not already in ring), add it to the grown set.
                if nbr.GetAtomicNum() == 8:
                    grown.add(nbr.GetIdx())
        grown_sets.append(grown)
    
    # Build a connectivity graph among these 4 units.
    # Two units are considered connected if their grown sets share at least one atom.
    n = len(grown_sets)
    connectivity = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if grown_sets[i].intersection(grown_sets[j]):
                connectivity[i].add(j)
                connectivity[j].add(i)
    
    # Check that all 4 candidate units are mutually connected.
    visited = set()
    def dfs(node):
        stack = [node]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            stack.extend(list(connectivity[curr] - visited))
    
    dfs(0)
    if len(visited) != n:
        return False, "Candidate carbohydrate rings are not connected via glycosidic bonds"
    
    return True, "Molecule contains exactly 4 connected monosaccharide (sugar) units"

# For testing purposes (script usage)
if __name__ == '__main__':
    # Example tetrasaccharide from literature:
    test_smiles = "OC[C@H]1O[C@@H](OC[C@H]2O[C@@H](OC[C@H]3O[C@@H](OC[C@H]4OC(O)[C@H](O)[C@@H](O)[C@H]4O)[C@H](O)[C@@H](O)[C@H]3O)[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
    result, reason = is_tetrasaccharide(test_smiles)
    print(result, reason)
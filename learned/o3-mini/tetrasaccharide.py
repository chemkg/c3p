"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition: An oligosaccharide comprising four monomeric monosaccharide units.
The revised classifier:
  - Parses the SMILES and requires one connected component.
  - Finds rings of size 5 (furanose) or 6 (pyranose) that have exactly one oxygen in the ring.
  - For each candidate ring, it checks that there are exocyclic O–H substituents (at least one)
    and that a number of ring atoms are chiral (≥2 for 5-membered and ≥3 for 6-membered rings).
  - Finally, it confirms that the candidates are connected (via bonds between ring atoms)
    so that all four rings form one unified oligosaccharide.
If these tests pass, the molecule is classified as a tetrasaccharide.
"""

from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determine if a molecule (given as a SMILES string) is a tetrasaccharide.
    It checks:
      - Valid SMILES and a single connected component.
      - Exactly 4 candidate sugar rings. A candidate sugar ring is defined as:
          * A 5-membered or 6-membered ring containing exactly one oxygen (and 4 or 5 carbons, respectively).
          * At least one exocyclic O–H substituent (i.e. an oxygen neighbour not in the ring that carries a hydrogen).
          * A minimum number of chiral centers in the ring (≥2 if 5-membered, ≥3 if 6-membered).
      - The 4 candidate sugar rings are connected via bonds (i.e. the sugar units are linked by glycosidic bonds).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): (True, reason) if molecule qualifies as a tetrasaccharide;
                     (False, reason) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check that molecule is a single connected component.
    mol_smiles = Chem.MolToSmiles(mol)
    if '.' in mol_smiles:
        return False, "Molecule contains multiple disconnected fragments"
    
    ring_info = mol.GetRingInfo().AtomRings()
    
    def is_sugar_ring(atom_indices):
        """
        Decide whether the ring given by atom_indices qualifies as a monosaccharide ring.
        Checks:
          - Ring size must be 5 or 6.
          - Exactly one oxygen in the ring.
          - Exocyclic –OH: at least one neighbor of a ring atom is an oxygen with at least one H.
          - Chiral centers in the ring: at least 2 for a 5-membered ring, or at least 3 for a 6-membered ring.
        """
        n_atoms = len(atom_indices)
        if n_atoms not in (5, 6):
            return False
        
        oxygen_in_ring = 0
        carbon_in_ring = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_in_ring += 1
            elif atom.GetAtomicNum() == 6:
                carbon_in_ring += 1
        # A typical pyranose ring is 6-membered (1 oxygen + 5 carbons),
        # and a typical furanose ring is 5-membered (1 oxygen + 4 carbons).
        if n_atoms == 6 and not (oxygen_in_ring == 1 and carbon_in_ring == 5):
            return False
        if n_atoms == 5 and not (oxygen_in_ring == 1 and carbon_in_ring == 4):
            return False

        # Check for exocyclic hydroxyl groups.
        # For each atom in the ring, look at neighbors not in the ring.
        exo_oh_count = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in atom_indices:
                    continue
                # Check if neighbor is oxygen and (explicitly or implicitly) bonded to at least one hydrogen.
                if nbr.GetAtomicNum() == 8:
                    # A simple test: if the neighbor has any explicit hydrogen or an atomic number hint from implicit H count.
                    if nbr.GetTotalNumHs() > 0:
                        exo_oh_count += 1
        if exo_oh_count < 1:
            return False

        # Check chiral centers in the ring.
        chiral_count = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            # If the atom has a specified chiral tag (different from CHI_UNSPECIFIED) we count it.
            if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                chiral_count += 1
        required_chiral = 2 if n_atoms == 5 else 3
        if chiral_count < required_chiral:
            return False

        return True

    # Gather indices for candidate sugar rings.
    candidate_rings = []
    for ring in ring_info:
        if is_sugar_ring(ring):
            candidate_rings.append(set(ring))
    
    if len(candidate_rings) != 4:
        return (False, f"Found {len(candidate_rings)} candidate carbohydrate rings, need exactly 4")
    
    # Check connectivity among candidate rings.
    # We use a simple graph where each node is one candidate ring.
    # An edge is added if any atom in one ring is bonded to any atom in the other.
    n = len(candidate_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # If any bond connects an atom in candidate_rings[i] and an atom in candidate_rings[j], mark them as connected.
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
    
    # Now do a simple connectivity search.
    visited = set()
    def dfs(node):
        stack = [node]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            stack.extend(graph[cur] - visited)
    
    dfs(0)
    if len(visited) != n:
        return (False, "Candidate carbohydrate rings are not connected via glycosidic bonds")
    
    return True, "Molecule contains exactly 4 connected monosaccharide (sugar) units"

# For testing (if running as script):
if __name__ == '__main__':
    # Test with a known tetrasaccharide SMILES
    test_smiles = "OC[C@H]1O[C@@H](OC[C@H]2O[C@@H](OC[C@H]3O[C@@H](OC[C@H]4OC(O)[C@H](O)[C@@H](O)[C@H]4O)[C@H](O)[C@@H](O)[C@H]3O)[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
    result, reason = is_tetrasaccharide(test_smiles)
    print(result, reason)
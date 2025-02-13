"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins (Natural pigments containing a fundamental skeleton of four pyrrole nuclei 
united through the alpha-positions by four methine groups to form a macrocyclic structure.)

The approach is to:
1. Parse the SMILES using RDKit.
2. Identify all aromatic five‐membered rings that contain exactly one nitrogen (i.e. pyrrole rings).
3. Build a “graph” where each pyrrole ring is a node and an edge exists if two rings share exactly one atom.
4. Look for a cycle of four fused pyrrole rings in that graph.
Note that this is a heuristic that works for many known porphyrin derivatives.
"""

from rdkit import Chem
import itertools

def is_porphyrins(smiles: str):
    """
    Determines if a molecule belongs to the porphyrin class based on its SMILES string.
    
    A porphyrin is characterized by four pyrrole rings (five-membered aromatic rings with one nitrogen)
    connected via meso (methine) carbons into a macrocyclic structure.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a porphyrin.
        str: Explanation of the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Find candidate pyrrole rings:
    pyrrole_rings = []
    for ring in ring_info:
        if len(ring) == 5:  # pyrrole rings are 5-membered
            # Check if all atoms in this ring are aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                # Count the number of nitrogen atoms in the ring.
                n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                if n_count == 1:
                    pyrrole_rings.append(ring)
    
    if len(pyrrole_rings) < 4:
        return False, f"Found only {len(pyrrole_rings)} aromatic five‐membered rings with one nitrogen (pyrrole candidates), need at least 4."
    
    # Build a graph among pyrrole rings:
    # Two rings are considered connected if they share exactly one atom (the bridging meso carbon).
    graph = {i: set() for i in range(len(pyrrole_rings))}
    for i in range(len(pyrrole_rings)):
        for j in range(i+1, len(pyrrole_rings)):
            shared_atoms = set(pyrrole_rings[i]).intersection(pyrrole_rings[j])
            if len(shared_atoms) == 1:
                graph[i].add(j)
                graph[j].add(i)
                
    # Look for a cycle of four fused pyrrole rings.
    # This means there is an ordering of four rings such that:
    #   ring A is connected to ring B, B to C, C to D, and D back to A.
    for combo in itertools.combinations(range(len(pyrrole_rings)), 4):
        # Check all permutations of these 4 nodes
        for perm in itertools.permutations(combo):
            if (perm[1] in graph[perm[0]] and
                perm[2] in graph[perm[1]] and
                perm[3] in graph[perm[2]] and
                perm[0] in graph[perm[3]]):
                return True, "Porphyrin macrocycle detected via 4 fused pyrrole rings."
    
    return False, "No closed cycle of 4 fused pyrrole rings (porphyrin core) detected."
    
# Example usage:
if __name__ == '__main__':
    # A known porphyrin-like SMILES (simplified porphine) could be tested here.
    test_smiles = "c1cc2ncccc2c2ccccn12"  # (This is a simplified representation and may not be an exact porphyrin)
    result, reason = is_porphyrins(test_smiles)
    print("Is porphyrin?", result)
    print("Reason:", reason)
"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI saccharolipid (Lipids that contain a carbohydrate moiety.)
Improved version:
  - For the lipid part we compute connected components of non‐ring carbon atoms 
    and require one component have at least 7 atoms.
  - For the sugar part we look for 5– or 6–membered rings composed only of C and O,
    that are not aromatic.
  - Finally we require that at least one sugar ring is directly connected (by at least one bond)
    to an atom within one of the long acyclic (lipid) components.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a molecule that contains a long lipid chain
    (here defined as a connected acyclic chain of at least 7 carbon atoms)
    that is covalently attached to a carbohydrate moiety (a non‐aromatic 5– or 6–membered ring
    composed solely of carbon and oxygen, with at least one oxygen present).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a saccharolipid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    #############################
    # Step 1. Identify candidate lipid chains.
    # We consider only carbon atoms that are not in any ring.
    non_ring_carbons = [atom.GetIdx() for atom in mol.GetAtoms() 
                          if atom.GetAtomicNum() == 6 and not atom.IsInRing()]

    # Build an undirected graph (as a dict) for these non-ring carbons:
    # graph[node] = set(neighbor nodes) (only if the neighbor is also in non_ring_carbons)
    graph = {idx: set() for idx in non_ring_carbons}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in graph: 
            continue
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in graph:
                graph[idx].add(nb_idx)
                graph[nb_idx].add(idx)
                
    # Get the connected components of the non-ring carbon graph.
    def connected_components(graph_dict):
        seen = set()
        comps = []
        for node in graph_dict:
            if node not in seen:
                # do depth-first search
                stack = [node]
                comp = set()
                while stack:
                    n = stack.pop()
                    if n not in seen:
                        seen.add(n)
                        comp.add(n)
                        stack.extend(graph_dict[n] - seen)
                comps.append(comp)
        return comps
    
    comps = connected_components(graph)
    # From all components, pick those that have at least 7 atoms.
    valid_lipid_components = [comp for comp in comps if len(comp) >= 7]
    if not valid_lipid_components:
        return False, "No long acyclic carbon chain (lipid chain) of at least 7 carbons detected"
    
    #############################
    # Step 2. Identify candidate sugar rings.
    # We inspect all rings (from GetRingInfo) and look for rings that:
    #   - Have size 5 or 6.
    #   - Contain only carbon and oxygen atoms.
    #   - Are not aromatic.
    ring_info = mol.GetRingInfo()
    sugar_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        # Get atomic numbers and aromaticity for atoms in the ring.
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        nums = [atom.GetAtomicNum() for atom in atoms]
        # Check that every atom is either C (6) or O (8) and at least one oxygen is present.
        if all(num in (6, 8) for num in nums) and any(num == 8 for num in nums):
            # Also require that none of the ring atoms is aromatic.
            if not any(atom.GetIsAromatic() for atom in atoms):
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety (5- or 6-membered, non-aromatic ring of C/O) detected"
    
    #############################
    # Step 3. Check connectivity between a sugar ring and a lipid chain.
    # We will require that at least one atom from a sugar ring has a bond to an atom
    # that is part of one of the valid lipid (non-ring carbon) components.
    connected = False
    # For quick lookup: merge all valid lipid chain atom indices into a set.
    lipid_atoms = set()
    for comp in valid_lipid_components:
        lipid_atoms.update(comp)
    
    # Now, for each sugar ring we check if any of its atoms is directly bonded to an atom in lipid_atoms.
    for ring in sugar_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in lipid_atoms:
                    connected = True
                    break
            if connected:
                break
        if connected:
            break
    if not connected:
        return False, "Sugar (carbohydrate) moiety not covalently attached to a lipid chain"

    return True, "Molecule contains a long acyclic lipid chain attached to a carbohydrate moiety"

# Example usage:
if __name__ == "__main__":
    # An example saccharolipid (taken from one of the provided examples)
    smiles_example = "CCCCCCCCCCCCC(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC"
    result, reason = is_saccharolipid(smiles_example)
    print("Is saccharolipid?", result)
    print("Reason:", reason)
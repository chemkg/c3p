"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid
Definition:
  Any triterpenoid that is highly oxygenated and has a prototypical structure either containing or derived from 
  a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton. Examples include limonin, azadirachtin, and related compounds.
  
Heuristic criteria in this implementation:
  - Molecular weight between 300 and 820 Da.
  - Carbon atom count between 24 and 40.
  - At least 4 oxygen atoms in the molecule.
  - If oxygen count is moderate (4 to 6), the molecule must include a typical furan ring (SMARTS "c1ccoc1").
  - The fused polycyclic core (assessed via ring connectivity of rings sharing at least 2 atoms) should have at least 3 rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string using heuristic criteria.
    
    Args:
       smiles (str): SMILES representation of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a limonoid, otherwise False.
       str: Explanation of the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (300 - 820 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 820:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the expected range (300-820 Da)"
    
    # Check the count of carbon atoms (typical limonoids have 24-40 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 24 or c_count > 40:
        return False, f"Carbon count ({c_count}) is outside the range expected for a triterpenoid skeleton (24-40)"
    
    # Check oxygen count: limonoids are highly oxygenated (at least 4 oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Not enough oxygen atoms ({o_count}); expected at least 4 for a limonoid"

    # Check for fused polycyclic core: get all ring atom lists from RDKit
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found; limonoids are polycyclic"
    
    # Build a graph where each ring is a node and add an edge if two rings share at least two atoms
    n_rings = len(ring_info)
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(set(ring_info[i]).intersection(ring_info[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
                
    # Find the largest connected component (fused ring network) in the ring graph
    visited = set()
    def dfs(node):
        stack = [node]
        comp = set()
        while stack:
            cur = stack.pop()
            if cur not in comp:
                comp.add(cur)
                for neighbor in adj[cur]:
                    if neighbor not in comp:
                        stack.append(neighbor)
        return comp

    max_fused = 0
    for i in range(n_rings):
        if i not in visited:
            comp = dfs(i)
            visited.update(comp)
            if len(comp) > max_fused:
                max_fused = len(comp)
                
    # Relaxed threshold: expect at least 3 fused rings for the limonoid core
    if max_fused < 3:
        return False, f"Fused ring network too small ({max_fused} fused rings found); expected at least 3 in a limonoid core"
    
    # Check for the presence of a furan ring if oxygen count is moderate (4 to 6) 
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    has_furan = mol.HasSubstructMatch(furan_pattern)
    if 4 <= o_count <= 6 and not has_furan:
        return False, "Furan ring not found (oxygenation is only moderate); expected a typical furan ring in limonoids"
    
    return True, ("Molecule appears to be a limonoid (triterpenoid core with 24-40 carbons, "
                  f"{o_count} oxygens, fused polycyclic core with {max_fused} rings, "
                  f"{'with' if has_furan else 'without'} a furan ring as needed)")
    
# Example usage (uncomment to test):
# test_smiles = "COC(=O)C[C@H]1[C@]2(C)C[C@@]3(O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](OC(=O)C(C)(C)O)C4=C1[C@H](OC(C)=O)[C@@]3(O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1"
# result, reason = is_limonoid(test_smiles)
# print(result, reason)
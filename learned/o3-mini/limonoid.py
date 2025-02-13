"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid
Definition:
  Any triterpenoid that is highly oxygenated and has a prototypical structure either containing or derived from 
  a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton. Examples include limonin, azadirachtin, and related compounds.
  
This revised implementation uses heuristic criteria:
  - Molecular weight must be between 300 and 820 Da.
  - Carbon count must be roughly between 24 and 40.
  - The molecule must have at least 4 oxygen atoms.
  - If oxygenation is moderate (4 to 6 oxygen atoms), a furan ring (SMARTS "c1ccoc1") must be present.
  - The molecule must have a fused polycyclic core: the largest connected (fused) ring network must contain at least 4 rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    
    Heuristic strategy:
      - Parse the SMILES.
      - Check that the molecular weight is between 300 and 820 Da.
      - Verify that the number of carbon atoms is roughly between 24 and 40.
      - Require that the molecule contains at least 4 oxygen atoms.
      - If oxygen count is moderate (4 to 6), require a furan ring substructure.
      - Ensure the molecule has a fused polycyclic core by requiring that the largest group of rings (that share
        at least 2 atoms with one another) contains at least 4 rings.
      
    Args:
       smiles (str): SMILES representation of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a limonoid, otherwise False.
       str: Explanation of the classification decision.
    """
    # Parse SMILES into an rdkit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Molecular weight criterion (relaxed upper bound: 300 - 820 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 820:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the expected range (300-820 Da)"
    
    # Carbon count criterion (limonoids typically have 24-40 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 24 or c_count > 40:
        return False, f"Carbon count ({c_count}) is outside the typical range for a triterpenoid skeleton (24-40)"
    
    # Oxygen count (must be at least 4)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Not enough oxygen atoms ({o_count}); expected at least 4 for a limonoid"
    
    # Check for polycyclic core by assessing fused ring network
    # Get all rings as tuples of atom indices.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found; limonoids are polycyclic"
    
    # Build an adjacency graph where each node is a ring (by its index)
    # Two rings are considered fused if they share at least two atoms
    n_rings = len(ring_info)
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            set_i = set(ring_info[i])
            set_j = set(ring_info[j])
            if len(set_i.intersection(set_j)) >= 2:
                adj[i].add(j)
                adj[j].add(i)
                
    # Find the size of the largest connected component in the ring graph.
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
            visited |= comp
            if len(comp) > max_fused:
                max_fused = len(comp)
    if max_fused < 4:
        return False, f"Fused ring network too small ({max_fused} fused rings found); expected at least 4 in a limonoid core"
    
    # Check for furan ring: if oxygenation is moderate (4 to 6 oxygen atoms)
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    has_furan = mol.HasSubstructMatch(furan_pattern)
    if 4 <= o_count <= 6 and not has_furan:
        return False, "Furan ring not found (and oxygenation is only moderate); expected typical furan ring in limonoids"
    
    return True, "Molecule appears to be a limonoid (triterpenoid core, highly oxygenated, fused polycyclic structure and, if oxygenation is moderate, contains a furan ring)"

# Example usage (uncomment to test):
# test_smiles = "COC(=O)C[C@H]1[C@]2(C)C[C@@]3(O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](OC(=O)C(C)(C)O)C4=C1[C@H](OC(C)=O)[C@@]3(O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1"
# result, reason = is_limonoid(test_smiles)
# print(result, reason)
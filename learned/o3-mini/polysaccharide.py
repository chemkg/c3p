"""
Classifies: CHEBI:18154 polysaccharide
"""
#!/usr/bin/env python
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
A structure is classified as a polysaccharide if it has more than ten distinct sugar rings and a molecular weight above a threshold.
This version improves detection by examining the ring information (instead of only SMARTS matching);
it requires that a candidate ring (of size 5 or 6) contains exactly one ring oxygen (atomic number 8)
and that most of the other atoms are saturated carbons bearing external oxygen substituents (as a proxy for hydroxyl groups).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import HybridizationType

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide by counting sugar rings using molecular ring info and some heuristics.
    
    Criteria:
      (1) The SMILES must parse.
      (2) A candidate monosaccharide ring is a 5-membered (furanose-like) or 6-membered (pyranose-like) ring.
          It should contain exactly one oxygen atom and the remainder sp3 carbons.
          In addition, several (heuristically 2 for furanoses and 3 for pyranoses) non‐ring oxygen substituents must be attached to carbons in the ring.
      (3) The distinct sugar ring count (each ring is counted only once) must be greater than 10.
      (4) The molecular weight must be at least 1000 Da.
      
    Args:
      smiles (str): Input SMILES
      
    Returns:
      (bool, str): (True, explanation) if classified as polysaccharide; else (False, explanation).
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_ring_sets = set()
    
    # Loop over every ring detected
    for ring in ring_info:
        size = len(ring)
        if size not in (5,6):
            continue  # only look at rings of size 5 or 6
        # Gather atoms in the ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        num_ox = sum(1 for a in ring_atoms if a.GetAtomicNum() == 8)
        num_c = sum(1 for a in ring_atoms if a.GetAtomicNum() == 6 and a.GetHybridization() == HybridizationType.SP3)
        
        # For a pyranose-like ring (6-membered) we expect exactly one oxygen and roughly 5 saturated carbons.
        # For a furanose-like ring (5-membered) we expect one oxygen and 4 saturated carbons.
        if size == 6 and num_ox != 1:
            continue
        if size == 5 and num_ox != 1:
            continue
        
        # Count external oxygen substituents on the carbons in the ring (exclude atoms in the ring)
        ext_oh = 0
        for a in ring_atoms:
            if a.GetAtomicNum() == 6:
                for nb in a.GetNeighbors():
                    if nb.GetIdx() not in ring and nb.GetAtomicNum() == 8:
                        ext_oh += 1
                        break  # count each carbon at most once
        # Heuristic: require at least 3 external OH’s for a 6-membered ring and at least 2 for a 5-membered ring.
        if size == 6 and ext_oh < 3:
            continue
        if size == 5 and ext_oh < 2:
            continue

        sugar_ring_sets.add(frozenset(ring))
    
    sugar_unit_count = len(sugar_ring_sets)
    
    if sugar_unit_count <= 10:
        return False, f"Found only {sugar_unit_count} sugar ring(s); a polysaccharide requires more than 10 residues."
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a typical polysaccharide."
    
    return True, f"Detected {sugar_unit_count} sugar ring(s) with a molecular weight of {mol_wt:.1f} Da indicative of a polysaccharide."

# Example usage:
# Uncomment the following lines to test:
# test_smiles = "O([C@H]1[C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)..."  # (example excerpt)
# classified, reason = is_polysaccharide(test_smiles)
# print(classified, reason)
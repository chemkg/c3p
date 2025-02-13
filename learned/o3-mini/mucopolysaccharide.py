"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
#!/usr/bin/env python
"""
Classifies: mucopolysaccharide (glycosaminoglycan)
Definition: any of the group of polysaccharides composed of alternating units from
uronic acids and glycosamines, and commonly partially esterified with sulfuric acid.

This version uses several heuristics:
  1. The molecule must parse from the SMILES string.
  2. It must have at least three rings – and at least two rings should have the size 
     (5–6 atoms) and high oxygen content typical for sugar rings.
  3. It must have at least one uronic acid feature. We look for a ring carbon that 
     bears an exocyclic carboxyl group. (SMARTS: "[C;R]-[C;!R](=O)[O;H,O-]")
  4. It must have at least one glycosamine feature – a ring carbon that is substituted 
     with a primary amine. (SMARTS: "[C;R]-[N;H2]")
  5. The counts of these two features should be nearly equal (difference ≤ 1) to be 
     consistent with an alternating pattern.
  6. There must be at least one sulfate group (SMARTS: "[S](=O)(=O)[O-]") as many mucopolysaccharides are sulfated.
  7. The overall oxygen fraction should be high (≥ 0.33) as oxygen is abundant in sugars.
If any test fails, the function returns False with a short reason.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan) based on its SMILES.
    
    Heuristics:
      1. The SMILES must be valid.
      2. The molecule must include at least three rings.
      3. It should include at least two rings of typical size (5–6 atoms) and having at least 2 oxygens,
         to act as a proxy for sugar rings.
      4. It should have at least one uronic acid feature:
            A ring carbon covalently linked to an exocyclic carboxyl group.
         (SMARTS: "[C;R]-[C;!R](=O)[O;H,O-]")
      5. It should have at least one glycosamine feature:
            A ring carbon linked to a primary amine.
         (SMARTS: "[C;R]-[N;H2]")
      6. These two features must be nearly equal in count (difference ≤ 1).
      7. The molecule must show at least one sulfate group: (SMARTS: "[S](=O)(=O)[O-]").
      8. The oxygen fraction should be high (≥ 0.33).
    
    Args:
      smiles (str): SMILES string representation of the molecule.
    
    Returns:
      bool: True if it qualifies as a mucopolysaccharide, False otherwise.
      str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a minimum number of rings.
    ring_info = mol.GetRingInfo()
    total_rings = ring_info.NumRings()
    if total_rings < 3:
        return False, f"Too few rings detected ({total_rings}); expected at least 3 for a sugar polymer"
    
    # 2. Count candidate sugar rings: rings of size 5 or 6 with at least 2 oxygen atoms.
    sugar_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) in (5,6):
            oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_in_ring >= 2:
                sugar_ring_count += 1
    if sugar_ring_count < 2:
        return False, f"Too few candidate sugar rings detected ({sugar_ring_count}); expected at least 2"
    
    # 3. Uronic acid feature: a ring carbon with an exocyclic carboxyl group.
    # SMARTS: a ring carbon ([C;R]) linked to a non‐ring carbon that is a carbonyl ([C;!R](=O)[O;H,O-])
    uronic_smarts = "[C;R]-[C;!R](=O)[O;H,O-]"
    uronic_pattern = Chem.MolFromSmarts(uronic_smarts)
    if uronic_pattern is None:
        return False, "Failed to compile uronic acid SMARTS pattern"
    uronic_matches = mol.GetSubstructMatches(uronic_pattern)
    count_uronic = len(uronic_matches)
    if count_uronic < 1:
        return False, "No uronic acid (exocyclic carboxyl on a ring carbon) feature detected"
    
    # 4. Glycosamine feature: a ring carbon linked to a primary amine.
    glyco_smarts = "[C;R]-[N;H2]"
    glyco_pattern = Chem.MolFromSmarts(glyco_smarts)
    if glyco_pattern is None:
        return False, "Failed to compile glycosamine SMARTS pattern"
    glyco_matches = mol.GetSubstructMatches(glyco_pattern)
    count_glyco = len(glyco_matches)
    if count_glyco < 1:
        return False, "No glycosamine (primary amine on a ring carbon) feature detected"
    
    # 5. Check that the counts of uronic and glycosamine features are nearly equal (difference ≤ 1)
    if abs(count_uronic - count_glyco) > 1:
        return False, f"Feature counts not nearly equal: {count_uronic} uronic vs {count_glyco} glycosamine"
    
    # 6. Must be partially sulfated – look for at least one sulfate group.
    sulfate_smarts = "[S](=O)(=O)[O-]"
    sulfate_pattern = Chem.MolFromSmarts(sulfate_smarts)
    if sulfate_pattern is None:
        return False, "Failed to compile sulfate SMARTS pattern"
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if len(sulfate_matches) < 1:
        return False, "No sulfate groups detected; mucopolysaccharides are commonly sulfated"
    
    # 7. Check the overall oxygen fraction.
    atoms = list(mol.GetAtoms())
    if not atoms:
        return False, "No atoms in molecule"
    oxygen_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    oxy_fraction = oxygen_count / len(atoms)
    if oxy_fraction < 0.33:
        return False, f"Oxygen fraction too low ({oxy_fraction:.2f}); expected ≥ 0.33 for sugar polymers"
    
    msg = (f"Detected {total_rings} rings (with {sugar_ring_count} candidate sugar rings), "
           f"{count_uronic} uronic acid features and {count_glyco} glycosamine features (nearly equal), "
           f"plus sulfate groups and high oxygen fraction ({oxy_fraction:.2f}).")
    return True, msg

# Example usage for testing:
if __name__ == "__main__":
    # One of the provided examples (e.g. Desferrioxamine X1, expected to be a mucopolysaccharide by our definition)
    examples = [
        "ON1CCCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC1=O",  # Desferrioxamine X1
        "O=C1N2[C@H](C=CC=CC=CC=CC([C@@](C=CC=CC=C1)(O)C)O[C@@H]3OC[C@@H](O)[C@@H]([C@H]3N)O)[C@H](O)[C@H](C2)C",  # Ciromicin A
    ]
    for smi in examples:
        result, reason = is_mucopolysaccharide(smi)
        print(result, reason)
"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
#!/usr/bin/env python
"""
Classifies: mucopolysaccharide (glycosaminoglycan)
Definition: any polysaccharide composed of alternating units from uronic acids and glycosamines,
and commonly partially esterified with sulfuric acid.

This function uses relaxed heuristics:
  1. The molecule must parse from the SMILES string.
  2. It must contain at least one ring – typical sugar units are cyclic.
  3. It must have at least one uronic acid feature. Here we search for a ring carbon
     carrying a carbonyl and an oxygen substituent, where the oxygen can be protonated or deprotonated.
     (SMARTS: "[C;R](=O)[O;H,O-]")
  4. It must have at least one glycosamine feature. We look for a ring carbon bound to a primary amine.
     (SMARTS: "[C;R]-[N;H2]")
  5. We check that the counts of these two features are nearly equal (difference ≤ 1) to support an alternating pattern.
  6. We also require that the overall oxygen fraction is high (≥ 0.30) typical of oxygenated sugars.
  
If any of these tests fails, we return False along with a reason.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan) based on its SMILES.
    Heuristics:
      1. The molecule must be valid.
      2. It must have at least one ring (a minimal indication of cyclic sugar units).
      3. It should have at least one uronic acid feature (a ring carbon bound to a carbonyl and an -OH or -O– group).
      4. It should have at least one glycosamine feature (a ring carbon attached to a primary amine).
      5. The counts of the uronic acid and glycosamine features should be nearly equal (difference ≤ 1),
         supporting an alternating repeating unit.
      6. The molecule should be highly oxygenated (oxygen fraction ≥ 0.30).
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is classified as a mucopolysaccharide, False otherwise.
        str: A message explaining the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that at least one ring exists.
    if not mol.GetRingInfo().NumRings():
        return False, "No rings detected; unlikely to be a sugar-based polymer"

    # Define SMARTS for uronic acid feature.
    # The pattern looks for a ring carbon with a carbonyl group bonded to either a protonated or deprotonated oxygen.
    uronic_smarts = "[C;R](=O)[O;H,O-]"
    uronic_pattern = Chem.MolFromSmarts(uronic_smarts)
    if uronic_pattern is None:
        return False, "Failed to compile uronic acid SMARTS pattern"
    uronic_matches = mol.GetSubstructMatches(uronic_pattern)
    count_uronic = len(uronic_matches)
    
    # Define SMARTS for glycosamine feature.
    # This pattern matches a ring carbon attached to a nitrogen that has two hydrogens (primary amine).
    glyco_smarts = "[C;R]-[N;H2]"
    glyco_pattern = Chem.MolFromSmarts(glyco_smarts)
    if glyco_pattern is None:
        return False, "Failed to compile glycosamine SMARTS pattern"
    glyco_matches = mol.GetSubstructMatches(glyco_pattern)
    count_glyco = len(glyco_matches)
    
    if count_uronic < 1:
        return False, "No uronic acid (carboxylated sugar) moiety detected in a ring"
    if count_glyco < 1:
        return False, "No glycosamine (amino sugar) moiety detected in a ring"
    
    # Check if the two counts are nearly equal (allowing a difference of at most 1)
    if abs(count_uronic - count_glyco) > 1:
        return False, f"Counts not nearly equal (uronic: {count_uronic}, glycosamine: {count_glyco}); pattern not alternating"
    
    # Calculate oxygen fraction: ratio of oxygen atoms to all atoms.
    atoms = mol.GetAtoms()
    total_atoms = len(atoms)
    oxygen_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    oxy_fraction = oxygen_count / total_atoms if total_atoms > 0 else 0
    if oxy_fraction < 0.30:
        return False, f"Oxygen fraction too low ({oxy_fraction:.2f}); mucopolysaccharides are highly oxygenated"
    
    return True, ("The molecule contains ring systems with both uronic acid and glycosamine features in nearly equal numbers, "
                  "and exhibits a high oxygen content consistent with mucopolysaccharides.")

# Example usage for testing:
if __name__ == "__main__":
    # Test example: Desferrioxamine X1 (provided in the examples)
    test_smiles = "ON1CCCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC1=O"
    result, message = is_mucopolysaccharide(test_smiles)
    print(result, message)
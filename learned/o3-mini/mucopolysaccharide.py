"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
#!/usr/bin/env python
"""
Classifies: mucopolysaccharide (glycosaminoglycan)
Definition: any of the group of polysaccharides composed of alternating units from
uronic acids and glycosamines, and commonly partially esterified with sulfuric acid.

This updated function uses relaxed heuristics:
  1. The molecule must parse from the SMILES string.
  2. It must contain at least one ring – common for cyclic sugar units.
  3. It must have at least one uronic acid feature. Instead of requiring the carbonyl to be inside the ring,
     we now search for a pattern where a ring carbon is substituted with an exocyclic carboxyl group.
     (SMARTS: "[C;R]C(=O)[O;H,O-]")
  4. It must have at least one glycosamine feature – a ring carbon bearing an exocyclic primary amine.
     (SMARTS: "[C;R][N;H2]")
  5. We check that the counts of these two features are nearly equal (difference ≤ 1) to support an alternating pattern.
  6. We also require that the overall oxygen fraction is high (≥ 0.30) typical of oxygenated sugars.
  
If any test fails, we return False along with a brief reason.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan) based on its SMILES.
    
    Heuristics:
      1. The SMILES must be valid.
      2. The molecule must include at least one ring (sugar rings are cyclic).
      3. It should have at least one uronic acid feature – detected as a ring carbon substituted with a carboxyl group.
         (SMARTS: "[C;R]C(=O)[O;H,O-]")
      4. It should have at least one glycosamine feature – a ring carbon substituted with a primary amine.
         (SMARTS: "[C;R][N;H2]")
      5. The counts of these two features should be nearly equal (difference ≤ 1) to support an alternating unit pattern.
      6. The overall oxygen fraction should be high (≥ 0.30), as sugars are oxygen rich.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
      bool: True if the molecule passes all tests as a mucopolysaccharide.
      str: A message explaining the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of at least one ring (required for sugar rings)
    if not mol.GetRingInfo().NumRings():
        return False, "No rings detected; unlikely to be a sugar-based polymer"
    
    # Improved SMARTS for uronic acid: a ring carbon substituted by a carboxyl group (-C(=O)[O;H,O-])
    uronic_smarts = "[C;R]C(=O)[O;H,O-]"
    uronic_pattern = Chem.MolFromSmarts(uronic_smarts)
    if uronic_pattern is None:
        return False, "Failed to compile uronic acid SMARTS pattern"
    uronic_matches = mol.GetSubstructMatches(uronic_pattern)
    count_uronic = len(uronic_matches)
    
    # Improved SMARTS for glycosamine: a ring carbon substituted by a primary amine (-[N;H2])
    glyco_smarts = "[C;R][N;H2]"
    glyco_pattern = Chem.MolFromSmarts(glyco_smarts)
    if glyco_pattern is None:
        return False, "Failed to compile glycosamine SMARTS pattern"
    glyco_matches = mol.GetSubstructMatches(glyco_pattern)
    count_glyco = len(glyco_matches)
    
    # Check for presence of at least one uronic acid feature
    if count_uronic < 1:
        return False, "No uronic acid (carboxylated sugar) feature detected (expected as an exocyclic carboxyl on a ring carbon)"
    
    # Check for presence of at least one glycosamine feature
    if count_glyco < 1:
        return False, "No glycosamine (amino sugar) feature detected (expected as a primary amine on a ring carbon)"
    
    # Check that the counts are nearly equal (difference ≤ 1) to support an alternating pattern
    if abs(count_uronic - count_glyco) > 1:
        return False, f"Counts not nearly equal (uronic: {count_uronic}, glycosamine: {count_glyco}); pattern not alternating"
    
    # Calculate the overall oxygen fraction to verify that the molecule is highly oxygenated
    atoms = mol.GetAtoms()
    total_atoms = len(atoms)
    oxygen_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    oxy_fraction = oxygen_count / total_atoms if total_atoms else 0
    if oxy_fraction < 0.30:
        return False, f"Oxygen fraction too low ({oxy_fraction:.2f}); mucopolysaccharides are expected to be highly oxygenated"
    
    return True, ("The molecule exhibits ring systems with both exocyclic uronic acid and glycosamine features in nearly equal numbers, "
                  "and shows a high oxygen content consistent with mucopolysaccharides.")

# Example usage for testing:
if __name__ == "__main__":
    # For demonstration, testing one of the provided SMILES examples (Desferrioxamine X1)
    test_smiles = "ON1CCCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC1=O"
    result, message = is_mucopolysaccharide(test_smiles)
    print(result, message)
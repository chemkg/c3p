"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
#!/usr/bin/env python
"""
Classifies: mucopolysaccharide (glycosaminoglycan)
Definition: any polysaccharide composed of alternating units from uronic acids and glycosamines,
and commonly partially esterified with sulfuric acid.
This heuristic function checks if a molecule (given by its SMILES) meets some basic features
of such biopolymers.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan) based on its SMILES.
    Heuristics:
      1. The molecule must be valid.
      2. It should have at least 2 rings (to indicate a polymeric/sugar-chain system).
      3. It should contain both uronic acid and glycosamine features.
         - Uronic acids are sugar rings carrying a carboxylic acid; we search for a ring carbon
           with a double-bonded oxygen and an –OH (SMARTS: "[C;R](=O)[O;H]").
         - Glycosamines are sugar rings with an amino substitution; we search for a ring carbon 
           attached to an NH2 group (SMARTS: "[C;R]-[NH2]").
      4. Because glycosaminoglycans are highly oxygenated, we check the oxygen fraction in the whole molecule.
         A rough threshold of 30% oxygen atoms is used.
      5. Finally, we require the counts of the two substructures to be nearly equal (difference ≤1)
         so that an alternating pattern is plausible.
      
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if the molecule is classified as a mucopolysaccharide, False otherwise.
        str: A reason message for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that there are at least 2 rings (indication of polymeric sugar units)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Not enough rings to indicate a polysaccharide structure"
    
    # Define SMARTS for uronic acid feature: a ring carbon with a carboxylic acid moiety.
    uronic_smarts = "[C;R](=O)[O;H]"
    uronic_pattern = Chem.MolFromSmarts(uronic_smarts)
    if uronic_pattern is None:
        return False, "Failed to compile uronic acid SMARTS pattern"
    uronic_matches = mol.GetSubstructMatches(uronic_pattern)
    count_uronic = len(uronic_matches)
    
    # Define SMARTS for glycosamine feature: a ring carbon attached to an –NH2 group.
    glyco_smarts = "[C;R]-[NH2]"
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
    
    # Calculate oxygen fraction in the molecule.
    atoms = mol.GetAtoms()
    total_atoms = len(atoms)
    oxygen_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    oxy_fraction = oxygen_count / total_atoms if total_atoms > 0 else 0
    if oxy_fraction < 0.30:
        return False, f"Oxygen fraction too low ({oxy_fraction:.2f}); mucopolysaccharides are highly oxygenated"
    
    return True, ("The molecule shows multiple rings with features consistent with alternating "
                  "uronic acid (carboxylated) and glycosamine (amino sugar) units and has high oxygen content.")

# Example usage:
if __name__ == "__main__":
    # Try one of the provided example SMILES (for testing purposes)
    test_smiles = "ON1CCCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC(=O)N(O)CCCCNC(=O)CCC1=O"  # Desferrioxamine X1
    result, message = is_mucopolysaccharide(test_smiles)
    print(result, message)
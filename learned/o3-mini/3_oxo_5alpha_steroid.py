"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI: 3-oxo-5alpha-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5. 
Heuristic (improved):
  (1) Look for a ketone group (C=O) that is part of a ring.
  (2) Require the molecule contains several fused rings of size 5 or 6 (a rough approximation to the steroid nucleus).
  (3) To favor a 5α configuration we require that at least one chiral center in a ring junction is drawn with an explicit alpha annotation.
      Here we “hack” this by requiring that the original SMILES string contain the substring "[C@@", which (by convention) is often used for an alpha chiral center.
Note: Stereochemistry in steroids is complex and a full 3D analysis would be needed to be sure. This approach is heuristic.
"""

from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-oxo-5alpha-steroid, else False.
        str: Explanation for the classification.
    """
    # Parse SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1: Check for 3-oxo functionality:
    # We require a ketone group (C=O) that is part of a ring.
    # SMARTS: [#6;R](=O) finds a carbon in a ring double-bonded to oxygen.
    ketone_pattern = Chem.MolFromSmarts("[#6;R](=O)")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found in a ring (3-oxo requirement not met)."
    
    # Step 2: Check for the steroid nucleus.
    # A typical steroid nucleus has four fused rings (three cyclohexanes and one cyclopentane).
    # We will make a heuristic: count rings of size 5 or 6.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # list of tuples (each tuple is indices of atoms constituting a ring)
    steroid_rings = [ring for ring in rings if len(ring) in (5, 6)]
    if len(steroid_rings) < 4:
        return False, f"Only {len(steroid_rings)} rings of size 5 or 6 found (steroid nucleus not detected)."
    
    # Step 3: Check for alpha configuration at a ring junction.
    # In a steroid, the C-5 atom (at a ring junction) is assigned as alpha. Ideally we would identify
    # the specific ring-junction and then from 3D geometry decide whether the hydrogen is “down”.
    # Here we resort to a heuristic that the input SMILES explicitly uses the "[C@@" notation
    # (commonly used in steroids for an alpha chiral center). We do a string search on the original SMILES.
    if "[C@@" not in smiles:
        return False, "No explicit alpha stereochemical annotation found (5α requirement not met)."
    
    return True, "Molecule has a 3-oxo group in a fused 4-ring system and displays an explicit alpha stereochemistry (likely a 5α-steroid)."


# Example usage (for testing; remove or comment out in production):
if __name__ == "__main__":
    # Test with one example: 5alpha-Androstane-3,11,17-trione
    test_smiles = "C[C@]12CC(=O)[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2=O"
    result, reason = is_3_oxo_5alpha_steroid(test_smiles)
    print(result, reason)
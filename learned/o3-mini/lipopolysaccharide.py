"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
#!/usr/bin/env python
"""
Classifies: A lipopolysaccharide (LPS)

Definition (heuristic): Natural compounds consisting of a trisaccharide (repeating unit)
composed of two heptose units and one octulosonic acid (a KDO sugar) plus oligosaccharide side‐chains
and one or more 3-hydroxytetradecanoic acid units (as found in Gram-negative bacterial LPS).
Due to the structural complexity, the approach is heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    
    Heuristic approach:
      - The molecule must be a valid structure.
      - Its molecular weight should be reasonably high (here we require > 500 Da).
      - It should contain several "sugar-like" rings. Here we define a sugar ring as a 
        5-membered or 6-membered ring containing exactly one oxygen (as in a typical furanose or pyranose).
        We require at least 3 such rings (which roughly corresponds to the trisaccharide repeating unit
        and additional sugar side-chains).
      - It must contain at least one fatty acid substructure corresponding roughly to 
        3-hydroxytetradecanoic acid. We use the SMARTS pattern "O=C(O)CCCC(O)CCCCCCCC" as an approximation.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is consistent with a lipopolysaccharide, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES and check validity
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (LPS are usually heavy molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low for a lipopolysaccharide (wt = {mol_wt:.1f} Da)"
    
    # Count sugar-like rings:
    # We assume that a sugar ring is a 5- or 6-membered ring containing exactly one oxygen atom.
    sugar_count = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            oxygens_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O")
            if oxygens_in_ring == 1:
                sugar_count += 1
    if sugar_count < 3:
        return False, f"Only {sugar_count} sugar-like ring(s) detected; at least 3 are expected in a lipopolysaccharide"
    
    # Search for a 3-hydroxytetradecanoic acid unit.
    # This is a heuristic SMARTS pattern.
    fatty_acid_smarts = "O=C(O)CCCC(O)CCCCCCCC"
    fatty_acid_pattern = Chem.MolFromSmarts(fatty_acid_smarts)
    if fatty_acid_pattern is None:
        return False, "Could not parse fatty acid SMARTS pattern (internal error)"
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No 3-hydroxytetradecanoic acid-like substructure detected"
    
    # If all tests pass, we consider the molecule consistent with a lipopolysaccharide.
    return True, "Molecule has sufficient sugar-like rings and a fatty acid unit consistent with lipopolysaccharide structure"

# For testing purposes - you can remove or comment these lines if using this as a module.
if __name__ == "__main__":
    # Example: one of the provided SMILES strings (you may replace with specific examples)
    test_smiles = "O1C(OC(CC)C(O)=O)C(O)C(O)C(O)C1C(O)=O"  # 6-(1-carboxypropoxy)-3,4,5-trihydroxyoxane-2-carboxylic acid
    result, reason = is_lipopolysaccharide(test_smiles)
    print(result, reason)
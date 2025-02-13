"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
#!/usr/bin/env python
"""
Classifies: A lipopolysaccharide (LPS)

Definition (heuristic): LPS are natural compounds composed of a trisaccharide repeating unit 
(two heptose units and one octulosonic acid) with oligosaccharide side‐chains and one or more 
3-hydroxytetradecanoic acid units (as found in Gram-negative bacteria). Because of the 
structural complexity, we use a heuristic that checks for:
  - A valid chemical structure.
  - A molecular weight above ~350 Da.
  - The presence of multiple “sugar‐like” rings (we count rings having 5–7 atoms with one oxygen).
  - At least one fatty acid-like fragment (we search for a pattern approximating 3‑hydroxytetradecanoic acid).
  
Note: This heuristic will not always be perfect.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    
    Heuristic approach:
      - The SMILES string must parse to a molecule.
      - The molecular weight must be over ~350 Da (some known examples are as light as 396 Da).
      - The structure must contain several sugar-like rings. Here we scan for rings
        of size 5, 6, or 7 that have exactly one oxygen atom (this is a rough approximation
        to a pyranose, furanose, or heptose ring). We require at least 2 such rings.
      - The molecule must contain at least one fatty acid-like substructure. Here we define
        a substructure roughly matching a 3-hydroxytetradecanoic acid unit as a carboxylic acid 
        (O=C(O)) that is attached to a long (13-carbon) aliphatic chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is preferentially consistent with a lipopolysaccharide, False otherwise.
        str: A reason string.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight – lowering threshold here to capture lighter LPS candidates
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight too low for a lipopolysaccharide (wt = {mol_wt:.1f} Da)"
    
    # Count sugar-like rings.
    # Here we count rings having 5, 6, or 7 atoms that contain exactly one oxygen.
    sugar_count = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6, 7):
            oxygens_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O")
            if oxygens_in_ring == 1:
                sugar_count += 1
    # We require at least 2 sugar-like rings as a rough indicator of a trisaccharide (and possible side-chains)
    if sugar_count < 2:
        return False, f"Only {sugar_count} sugar-like ring(s) detected; at least 2 expected based on trisaccharide unit"
    
    # Search for a fatty acid-like fragment.
    # The following SMARTS pattern is a heuristic representation of a carboxylic acid attached to a long aliphatic chain.
    fatty_acid_smarts = "O=C(O)[CH2]CCCCCCCCCCCC"
    fa_pattern = Chem.MolFromSmarts(fatty_acid_smarts)
    if fa_pattern is None:
        return False, "Could not parse fatty acid SMARTS pattern (internal error)"
    if not mol.HasSubstructMatch(fa_pattern):
        return False, "No 3-hydroxytetradecanoic acid-like substructure detected"
        
    # If all tests pass, we consider the molecule consistent with a lipopolysaccharide.
    return True, "Molecule has features consistent with a lipopolysaccharide structure"

# For testing purposes – you may remove or comment out the code below if using this as a module.
if __name__ == "__main__":
    # Example test SMILES (you can replace with one of the provided examples)
    test_smiles = "O(C1[C@@H](O)[C@H](O)C(O[C@H]1OCCCC)COC(=O)CCCCCCCCCCCCCCC)[C@@H]2OC([C@H](OC(=O)CCC)[C@H](O)C2O)C"
    result, reason = is_lipopolysaccharide(test_smiles)
    print(result, reason)
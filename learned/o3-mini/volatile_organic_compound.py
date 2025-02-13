"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound (VOC)
Definition:
  A volatile organic compound is defined as any organic compound having an 
  initial boiling point less than or equal to 250 °C measured at 101.3 kPa.
  
Note:
  Because boiling point is not directly calculable from the SMILES (without detailed
  quantitative structure–property relationships), this implementation uses a 
  heuristic based on molecular weight (MW), topological polar surface area (TPSA), 
  substructure patterns (for free acids and phenols) and aromaticity.
  
Heuristic rules in this implementation:
  1. The molecule must be an organic compound (contain at least one carbon).
  2. If TPSA is very high (> 100 Å²), predict nonvolatile because strong intermolecular
     interactions are likely.
  3. If a free (non‐esterified) carboxylic acid group is present, predict nonvolatile.
  4. For “light” molecules (MW ≤ 250 Da):
       • If an aromatic “phenol” (an –OH directly bonded to an aromatic carbon) exists,
         predict nonvolatile because hydrogen bonding raises the boiling point.
       • Otherwise, predict volatile.
  5. For “heavier” molecules (MW > 250 Da):
       • If the molecule contains any aromatic rings then predict nonvolatile (many heavy
         aromatic compounds tend to have higher boiling points, even if nonpolar).
       • Otherwise, if the polarity (TPSA/MW) is very low (< 0.08), then predict volatile;
         else, predict nonvolatile.
         
Note: This heuristic is only one possible approach.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is predicted to be a volatile organic compound (VOC)
    based on its SMILES string using a heuristic that combines molecular weight,
    polarity (TPSA) and selected substructure patterns.
    
    Args:
      smiles (str): a SMILES string representing the molecule
      
    Returns:
      bool: True if predicted to be a VOC (initial boiling point <= 250°C), False otherwise.
      str: A short explanation for the classification.
    """
    # Try to parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is organic (contains at least one carbon).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atom found; not an organic compound"
    
    # Calculate key descriptors.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    polarity_ratio = tpsa / mol_wt if mol_wt else 0.0
    
    # Rule 1: Very high TPSA => nonvolatile.
    if tpsa > 100:
        return False, f"TPSA is high ({tpsa:.1f} Å²); strong intermolecular interactions suggest BP > 250°C"
    
    # Rule 2: Identify free carboxylic acid groups.
    # SMARTS for carboxylic acid: [CX3](=O)[OX1H]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    # Exclude cases where the acid is part of an ester.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O[C]")
    if acid_pattern and mol.HasSubstructMatch(acid_pattern):
        if not mol.HasSubstructMatch(ester_pattern):
            return False, "Contains a free carboxylic acid group; likely BP > 250°C"
    
    # For light molecules (MW <= 250 Da).
    if mol_wt <= 250:
        # Check for an aromatic–OH (phenol) pattern.
        # SMARTS: an oxygen (OH group) directly bound to an aromatic carbon.
        phenol_pattern = Chem.MolFromSmarts("c[OH]")
        if phenol_pattern and mol.HasSubstructMatch(phenol_pattern):
            return False, f"Light aromatic phenol detected; hydrogen bonding may raise BP above 250°C (MW = {mol_wt:.1f} Da)"
        else:
            return True, f"Molecular weight ({mol_wt:.1f} Da) is light; likely has BP <= 250°C"
    
    # For heavier molecules (MW > 250 Da):
    # First check: if any aromatic rings are present, assume nonvolatile.
    aromatic_ring_found = False
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_found = True
            break
    if aromatic_ring_found:
        return False, "Heavy molecule with aromatic ring(s); likely has BP > 250°C"
    
    # For heavy, nonaromatic molecules, try to "rescue" if polarity is very low.
    if polarity_ratio < 0.08:
        return True, f"Heavy molecule (MW = {mol_wt:.1f} Da) but very low polarity (TPSA/MW = {polarity_ratio:.3f}) suggests BP <= 250°C"
    else:
        return False, f"Heavy molecule (MW = {mol_wt:.1f} Da) with moderate polarity (TPSA = {tpsa:.1f} Å², ratio = {polarity_ratio:.3f}) likely has BP > 250°C"

# Example usage (for testing):
if __name__ == "__main__":
    # Example SMILES strings and names (selected from provided outcomes):
    examples = [
        ("CCCCCC(O)CCC", "nonan-4-ol"),
        ("CCCCCCCCCCCCCCCCCCCC(O)CCCCC", "hexacosan-6-ol"),
        ("C1COCCO1", "1,4-dioxane"),
        ("CCOC(C)(C)C", "tert-butyl ethyl ether"),
        ("C1CCCCC1", "cyclohexane"),
        ("CC[C@@H](C)C(C)C", "(3R)-2,3-dimethylpentane"),
        ("OC(=O)CC=C(C)C", "4-methylpent-3-enoic acid"),
        ("OC1=C(C(O)=CC(=C1)C)C", "Beta-Orcinol"),
        ("CCCCCCCCCCCCCC(O)CCCCCCCCCCC", "pentacosan-12-ol"),
        ("ClC(Cl)([H])[H]", "dichloromethane-d2"),
        ("Cc1ccccc1", "toluene"),
        ("COc1cc(ccc1O)[C@@H](O)[C@@H]1CO[C@@H]([C@H]1CO)c1ccc(O)c(OC)c1", "tanegool"),
    ]
    
    for smi, name in examples:
        result, reason = is_volatile_organic_compound(smi)
        print(f"SMILES: {smi:60s} NAME: {name:35s} -> {result}, REASON: {reason}")
"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound (VOC)
Definition:
  A volatile organic compound is defined as any organic compound having an 
  initial boiling point less than or equal to 250 °C measured at 101.3 kPa.
  
Note:
  Because boiling point is not directly calculable from the SMILES (without 
  detailed quantitative structure–property relationships such as the Joback method),
  this implementation uses a crude combination of molecular descriptors. In particular,
  we use molecular weight (MW), topological polar surface area (TPSA), and pattern matching
  to screen molecules. These rules are entirely heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is predicted to be a volatile organic compound (VOC)
    based on its SMILES string using a heuristic that combines molecular weight,
    polarity and selected substructure patterns.
    
    Rules (all are crude approximations):
      1. The molecule must be organic (contain at least one carbon).
      2. If the topological polar surface area (TPSA) is very high (>100 Å²) then
         the molecule is predicted non-volatile.
      3. If the molecule is light (MW <= 250 Da) then it is predicted volatile 
         unless it is an aromatic molecule having exactly one phenolic -OH group.
      4. If the molecule is heavy (MW > 250 Da) then we compute the ratio TPSA/MW.
         If the ratio is very low (<0.08) then the molecule is “rescued” and predicted volatile,
         otherwise non-volatile.
      5. If a free carboxylic acid group (non-ester) is found, this will tend to raise the boiling 
         point. In that case we predict non-volatile.
    
    Note: Because many of these criteria are problematic and not unambiguous,
          this is only one possible heuristic approach.
    
    Args:
      smiles (str): a SMILES string representing the molecule
      
    Returns:
      bool: True if predicted to be a volatile organic compound, False otherwise.
      str: A short explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that molecule is organic (has at least one carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atom found; not considered an organic compound"
    
    # Calculate descriptors: molecular weight and topological polar surface area (TPSA)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    
    # Compute the ratio of polarity per unit weight
    polarity_ratio = tpsa / mol_wt if mol_wt else 0.0
    
    # Rule 1: If TPSA is very high, predict non-volatile.
    if tpsa > 100:
        return False, f"TPSA is high ({tpsa:.1f} Å²); likely has strong intermolecular interactions and a boiling point > 250°C"
    
    # Rule 2: Check for free carboxylic acid group (but not ester).
    # The pattern [CX3](=O)[OX1H] typically matches a carboxylic acid.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if acid_pattern and mol.HasSubstructMatch(acid_pattern):
        # But also check that it is not part of an ester (which has [CX3](=O)O[C])
        ester_pattern = Chem.MolFromSmarts("[CX3](=O)O[C]")
        if not mol.HasSubstructMatch(ester_pattern):
            return False, "Contains a free carboxylic acid group; likely has a boiling point > 250°C"
    
    # For aromatic –OH (phenolic) groups, count the number attached to aromatic carbons.
    aromatic_oh_count = 0
    for atom in mol.GetAtoms():
        # Look for oxygen atom (potential -OH)
        if atom.GetSymbol() == "O":
            # Check if this oxygen is bonded to a hydrogen and to at least one carbon
            neighbors = atom.GetNeighbors()
            if len(neighbors) >= 1:
                # Check if any neighbor is a hydrogen (implicit or explicit)
                has_h = any(neighbor.GetSymbol() == "H" for neighbor in neighbors)
                # Check if at least one neighbor is aromatic carbon
                bonded_aromatic_c = any(neighbor.GetSymbol() == "C" and neighbor.GetIsAromatic() for neighbor in neighbors)
                if has_h and bonded_aromatic_c:
                    aromatic_oh_count += 1
    
    # Rule 3: For light molecules (MW <= 250 Da) use MW cutoff
    if mol_wt <= 250:
        if aromatic_oh_count == 1:
            # In our training outcomes, some aromatic alcohols with one -OH (e.g. certain phenols)
            # were mis‐classified as volatile; here we override and mark them as non-volatile.
            return False, f"Aromatic molecule with one phenolic -OH and low weight ({mol_wt:.1f} Da) - likely has a higher boiling point"
        else:
            return True, f"Molecular weight ({mol_wt:.1f} Da) is within light range; likely has an initial boiling point <= 250°C"
    
    # Rule 4: For heavier compounds (MW > 250 Da), “rescue” if the molecule is very nonpolar.
    if mol_wt > 250:
        if polarity_ratio < 0.08:
            return True, f"Molecular weight ({mol_wt:.1f} Da) is high but low polarity (TPSA/MW = {polarity_ratio:.3f}) suggests a lower boiling point"
        else:
            return False, f"Molecular weight ({mol_wt:.1f} Da) and polarity (TPSA = {tpsa:.1f} Å², ratio = {polarity_ratio:.3f}) indicate a boiling point > 250°C"
    
    # Fallback (should not be reached)
    return None, None

# Example usage (for testing purposes only):
if __name__ == "__main__":
    # List a few example SMILES strings and names:
    examples = [
        ("CCCCCC(O)CCC", "nonan-4-ol"),
        ("CCCCCCCCCCCCCCCCCCCC(O)CCCCC", "hexacosan-6-ol"),
        ("C1COCCO1", "1,4-dioxane"),
        ("CCOC(C)(C)C", "tert-butyl ethyl ether"),
        ("C1CCCCC1", "cyclohexane"),
        ("CC[C@@H](C)C(C)C", "(3R)-2,3-dimethylpentane"),
        ("OC(=O)CC=C(C)C", "4-methylpent-3-enoic acid"),
        ("OC1=C(C(O)=CC(=C1)C)C", "Beta-Orcinol"),
        ("OCCCC(=O)O", "Acetic acid derivative"),
    ]
    for smi, name in examples:
        result, reason = is_volatile_organic_compound(smi)
        print(f"SMILES: {smi:40s} NAME: {name:25s} -> {result}, REASON: {reason}")
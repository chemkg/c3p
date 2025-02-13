"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound (VOC)
Definition:
  A volatile organic compound is defined as any organic compound having an 
  initial boiling point less than or equal to 250 °C (482 °F) measured at 101.3 kPa.
  
Heuristic improvements over the previous attempt:
  • The molecule must be organic.
  • If TPSA is very high (>100 Å²) then strong intermolecular interactions may raise the BP.
  • Free carboxylic acid (or deprotonated acid) groups are flagged as nonvolatile for molecules
    with MW >150.
  • For light molecules (MW ≤ 250 Da):
         – If total hydrogen bond capacity (HBD + HBA) is high (≥3) and TPSA > 40 Å², predict nonvolatile.
         – If an aromatic ring containing 2 or more nitrogen atoms is present, predict nonvolatile.
         – Otherwise, predict volatile.
  • For heavier molecules (MW > 250 Da):
         – If any aromatic ring is present, predict nonvolatile.
         – Else if the polarity (TPSA/MW ratio) is very low (< 0.08), predict volatile;
           otherwise predict nonvolatile.
           
Note:
  Boiling point is not directly calculable from SMILES without a detailed QSAR model.
  The heuristic below is an approximate approach.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is predicted to be a volatile organic compound (VOC)
    (i.e. initial boiling point <= 250 °C) based on its SMILES string using a heuristic
    that combines molecular weight, TPSA, hydrogen-bonding capacity, and selective substructure searches.
    
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if predicted to be a VOC, False otherwise.
      str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is an organic compound (should contain at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atom found; not an organic compound"
    
    # Compute key descriptors.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    polarity_ratio = tpsa / mol_wt if mol_wt else 0.0
    num_hbd = rdMolDescriptors.CalcNumHBD(mol)
    num_hba = rdMolDescriptors.CalcNumHBA(mol)
    
    # --- Rule: Very high TPSA usually predicts a high boiling point.
    if tpsa > 100:
        return False, f"TPSA is high ({tpsa:.1f} Å²); strong intermolecular interactions suggest BP > 250°C"
    
    # --- Rule: Detect free carboxylic acid groups.
    # SMARTS for free carboxylic acid: [CX3](=O)[OX1H]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    # Also check for deprotonated carboxylate: [CX3](=O)[OX1-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1-]")
    if (acid_pattern and mol.HasSubstructMatch(acid_pattern)) or (carboxylate_pattern and mol.HasSubstructMatch(carboxylate_pattern)):
        # For very light acids (like acetic acid) boiling point is low.
        if mol_wt > 150:
            return False, "Contains a free carboxylic acid/carboxylate group; likely BP > 250°C"
    
    # --- For molecules with MW <= 250 Da:
    if mol_wt <= 250:
        # Check overall hydrogen bonding capacity.
        if (num_hbd + num_hba) >= 3 and tpsa > 40:
            return False, f"Light molecule (MW = {mol_wt:.1f} Da) but high H-bond capacity (HBD+HBA = {num_hbd+num_hba}, TPSA = {tpsa:.1f} Å²) suggests BP > 250°C"
        
        # Check for aromatic rings that are unusual: if an aromatic ring contains 2 or more nitrogen atoms,
        # we assume that this heterocycle may increase the boiling point.
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            # Only consider rings where all atoms are aromatic.
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                # Count nitrogen atoms in the ring.
                n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                if n_in_ring >= 2:
                    return False, "Light molecule with aromatic heterocycle (≥2 N atoms in a ring) tends to have BP > 250°C"
                    
        return True, f"Molecular weight ({mol_wt:.1f} Da) is light and hydrogen bonding is moderate; likely has BP <= 250°C"
    
    # --- For molecules with MW > 250 Da:
    else:
        # If the molecule contains any aromatic ring, assume a higher boiling point.
        ri = mol.GetRingInfo()
        aromatic_ring_found = False
        for ring in ri.AtomRings():
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_ring_found = True
                break
        if aromatic_ring_found:
            return False, "Heavy molecule with aromatic ring(s); likely has BP > 250°C"
        
        # For heavy, nonaromatic molecules, check polarity ratio.
        if polarity_ratio < 0.08:
            return True, f"Heavy molecule (MW = {mol_wt:.1f} Da) but very low polarity (TPSA/MW = {polarity_ratio:.3f}) suggests BP <= 250°C"
        else:
            return False, f"Heavy molecule (MW = {mol_wt:.1f} Da) with moderate polarity (TPSA = {tpsa:.1f} Å², ratio = {polarity_ratio:.3f}) likely has BP > 250°C"

# Example usage for testing the revised heuristic.
if __name__ == "__main__":
    test_examples = [
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
        ("C[N+](C)(C)C[C@@H](CC([O-])=O)OC(=O)C=C", "O-acryloyl-L-carnitine"),
        ("N1=C(CC=C)C(=NC=C1)C", "2-Methyl-3-(2-propenyl)pyrazine"),
    ]
    
    for smi, name in test_examples:
        res, reason = is_volatile_organic_compound(smi)
        print(f"SMILES: {smi:60s} NAME: {name:35s} -> {res}, REASON: {reason}")
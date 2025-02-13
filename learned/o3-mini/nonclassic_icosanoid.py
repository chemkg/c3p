"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: nonclassic icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
other than the classic icosanoids (the leukotrienes and the prostanoids).

Improved Heuristic rules:
  1. Molecule must be valid.
  2. Molecule must contain only allowed elements: C, H, O, N, S.
  3. Total number of carbon atoms between 18 and 26.
  4. Molecular weight between 250–500 Da.
  5. Must have at least 3 oxygen atoms and an oxygen/carbon ratio ≥ 0.3.
  6. Contains at least one carboxyl functionality (free acid, carboxylate, or ester).
  7. Contains at least three carbon–carbon double bonds (excluding C=O bonds).
  8. Excludes molecules with any 5‐membered ring.
  
Note: This heuristic approach is approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    
    Heuristic criteria:
      - Molecule is valid.
      - Allowed elements: C, H, O, N, and S.
      - Carbon count between 18 and 26.
      - Must have at least 3 oxygen atoms and oxygen/carbon ratio of at least 0.3.
      - Molecular weight between 250 and 500 Da.
      - Contains at least one carboxyl group (free acid, carboxylate, or ester).
      - Contains at least three C=C bonds between carbons (non-carbonyl double bonds).
      - Does not contain any 5-membered rings.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a nonclassic icosanoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Allow elements: H (1), C (6), N (7), O (8), S (16)
    allowed_atomic_nums = {1, 6, 7, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains element with atomic number {atom.GetAtomicNum()} not allowed. Allowed: C, H, O, N, S."
    
    # Count carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbons)
    if num_carbons < 18 or num_carbons > 26:
        return False, f"Carbon count {num_carbons} out of expected range (18–26) for a C20 precursor"
    
    # Count oxygen atoms.
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    num_oxygens = len(oxygens)
    if num_oxygens < 3:
        return False, f"Not enough oxygen atoms (found {num_oxygens}; need at least 3)"
    
    # Check oxygen/carbon ratio.
    if (num_oxygens/num_carbons) < 0.3:
        return False, f"Oxygen/carbon ratio ({num_oxygens}/{num_carbons} = {(num_oxygens/num_carbons):.2f}) below threshold (0.3)"
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da out of expected range (250–500 Da) for nonclassic icosanoids"
    
    # Check for carboxyl functionality.
    # Use SMARTS: free acid, carboxylate, and ester (as a proxy).
    carboxyl_smarts = ["C(=O)[OH]", "C(=O)[O-]", "C(=O)O[C]"]
    has_carboxyl = False
    for smarts in carboxyl_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            has_carboxyl = True
            break
    if not has_carboxyl:
        return False, "No carboxyl (acid/carboxylate/ester) group detected"
    
    # Count carbon–carbon double bonds (exclude carbonyl C=O bonds).
    alkene_count = 0
    for bond in mol.GetBonds():
        # Check if bond is double.
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Exclude if one of the atoms is oxygen (i.e. carbonyl).
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                alkene_count += 1
    if alkene_count < 3:
        return False, f"Not enough carbon–carbon double bonds (found {alkene_count}; need at least 3)"
    
    # Exclude molecules with any 5-membered rings.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            return False, "Contains a 5-membered ring – structure may be a classic prostanoid"
    
    return True, "Molecule meets heuristic criteria for nonclassic icosanoid"

# Optional testing code (remove or modify as needed).
if __name__ == "__main__":
    # Test examples (true positives and one false positive example provided)
    test_smiles = [
        "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O",  # true positive
        "CCCC\\C=C/C[C@@H]1O[C@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O",       # true positive
        "CCCC\\C=C/C\\C=C/CC(O)C(O)C\\C=C/CCCC([O-])=O"                    # false positive example
    ]
    
    for s in test_smiles:
        result, reason = is_nonclassic_icosanoid(s)
        print(result, reason)
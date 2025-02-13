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
  3. Total number of carbon atoms is expected to be close to 20 (here between 18 and 22).
  4. Molecular weight between 250–500 Da.
  5. Must have at least 3 oxygen atoms and an oxygen/carbon ratio of at least 0.20.
  6. Contains at least one carboxyl functionality (free acid or carboxylate).
  7. Contains at least three non-carbonyl C=C bonds.
  8. Excludes molecules with any aromatic atoms.
  9. Excludes molecules with any 5‐membered rings.
  
Note: This heuristic approach is approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    Heuristic criteria (improved):
      - Valid molecule.
      - Only allowed elements: C, H, O, N, and S.
      - Carbon count between 18 and 22.
      - Molecular weight between 250 and 500 Da.
      - At least 3 oxygen atoms and oxygen/carbon ratio >= 0.20.
      - Contains at least one carboxyl group (free acid or carboxylate).
      - Contains at least three carbon–carbon double bonds (excluding carbonyl bonds).
      - Contains no aromatic atoms.
      - Contains no 5-membered rings.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if classified as a nonclassic icosanoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Allowed elements: H (1), C (6), N (7), O (8), S (16)
    allowed_atomic_nums = {1, 6, 7, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains element with atomic number {atom.GetAtomicNum()} not allowed. Allowed: C, H, O, N, S."
    
    # Count carbon atoms and check range (close to 20)
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbons)
    if num_carbons < 18 or num_carbons > 22:
        return False, f"Carbon count {num_carbons} out of expected range (18–22) for a C20 precursor"
    
    # Count oxygen atoms.
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    num_oxygens = len(oxygens)
    if num_oxygens < 3:
        return False, f"Not enough oxygen atoms (found {num_oxygens}; need at least 3)"
    
    # Check oxygen/carbon ratio threshold (lowered to 0.2)
    ratio = num_oxygens/num_carbons
    if ratio < 0.20:
        return False, f"Oxygen/carbon ratio ({num_oxygens}/{num_carbons} = {ratio:.2f}) below threshold (0.20)"
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da out of expected range (250–500 Da) for nonclassic icosanoids"
    
    # Check for carboxyl functionality.
    # SMARTS to detect a free acid or carboxylate group.
    carboxyl_smarts = "[$([CX3](=O)[O;H,-])]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    if not (carboxyl_pattern and mol.HasSubstructMatch(carboxyl_pattern)):
        return False, "No carboxyl (acid/carboxylate) group detected"

    # Count non-carbonyl carbon–carbon double bonds.
    # For each double bond between carbon atoms, if neither carbon is in a carbonyl (i.e. double-bonded to oxygen) then count it.
    def is_carbonyl(atom):
        # returns True if the carbon atom has any double bond to oxygen.
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:
                    return True
        return False

    alkene_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only consider C=C bonds.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Exclude if either carbon is in a carbonyl.
                if not (is_carbonyl(a1) or is_carbonyl(a2)):
                    alkene_count += 1
    if alkene_count < 3:
        return False, f"Not enough non-carbonyl C=C bonds (found {alkene_count}; need at least 3)"
    
    # Exclude molecules with any aromatic atoms.
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            return False, "Contains aromatic atoms – not typical for a nonclassic icosanoid"
    
    # Exclude molecules with any 5-membered rings.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            return False, "Contains a 5-membered ring – structure may be a classic prostanoid or other artifact"
    
    return True, "Molecule meets heuristic criteria for nonclassic icosanoid"

# Optional testing code (remove or modify as needed).
if __name__ == "__main__":
    # Example test compounds (some true positives and known false positives are provided)
    test_smiles = {
        "20-hydroxylipoxin A4": "C(/C=C\\C=C\\[C@H](CCCCCO)O)=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)O",
        "20-hydroxylipoxin B4": "C(\\C=C/C=C/C=C/[C@H]([C@H](CCCCCO)O)O)=C/[C@H](CCCC(O)=O)O",
        # False positive examples – these should be rejected by our improvements.
        "Necroxime C": "O=C1O[C@H](C[C@@H](O)[C@H]2O[C@H]2[C@H](C=CC=3C1=C(O)C=CC3)O)CC=CNC(=O)/C=C/C=NOC",
        "20-carboxyleukotriene E4": "C(/C=C\\CCCCC(=O)O)/C=C\\C=C/C=C/[C@H]([C@@H](O)CCCC(=O)O)SC[C@H](N)C(=O)O"
    }
    
    for name, s in test_smiles.items():
        result, reason = is_nonclassic_icosanoid(s)
        print(f"{name}: {result} -> {reason}")
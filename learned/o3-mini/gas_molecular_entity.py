"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: gas molecular entity as defined by
"Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa)".

Heuristic changes from our previous approach:
  - Molecule must parse and be a single fragment.
  - If it is a single–atom species, only noble gases (He, Ne, Ar, Kr, Xe, Rn) are accepted.
  - All atoms must be from allowed main group elements, but we now exclude Si and B (which were leading to false positives).
  - Molecule must be small (molecular weight <=300 Da).
  - Molecule must be “rigid” (<=2 rotatable bonds).
  - The TPSA must be low; for molecules with weight <150 Da we allow TPSA < 60 Å², else TPSA must be < 40 Å².
  - Finally, if the molecule contains any rings then only 3-membered rings are allowed.
  
This heuristic approach may be imperfect.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string)
    qualifies as a gas molecular entity at STP.
    
    For our purposes, a gas is expected to:
      - Parse as a single–fragment molecule.
      - (If multi–atomic) contain only allowed main–group elements—with the proviso that
        compounds containing silicon or boron are excluded.
      - Have a relatively small molecular weight (<=300 Da).
      - Have low polar surface area: for molecules <150 Da, TPSA must be < 60 Å² and 
        for heavier molecules TPSA must be < 40 Å².
      - Be rigid (<=2 rotatable bonds).
      - Be structurally simple: if rings are present, only 3-membered rings (e.g. oxirane) are allowed.
      
    Single–atom molecules are only allowed if the atom is a noble gas (He, Ne, Ar, Kr, Xe, Rn).
    
    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
         bool: True if the molecule satisfies our criteria to be a gas, False otherwise.
         str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure we have a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (not a single molecular entity)"
    
    # Special handling for single–atom molecules.
    if mol.GetNumAtoms() == 1:
        # List of noble gas atomic numbers: He (2), Ne (10), Ar (18), Kr (36), Xe (54), Rn (86)
        noble_gases = {2, 10, 18, 36, 54, 86}
        atom = mol.GetAtomWithIdx(0)
        if atom.GetAtomicNum() in noble_gases:
            return True, "Single noble gas atom: gas at STP"
        else:
            return False, f"Single atom of element with atomic number {atom.GetAtomicNum()} is not a gas at STP"
    
    # Define allowed atomic numbers for main–group elements (from H up to Rn, excluding Si (14) and B (5)).
    allowed_atomic_numbers = {
         1,   2,    3, 4,   6, 7, 8, 9, 10,   # H, He, Li, Be, C, N, O, F, Ne
        11, 12,  13,    15, 16, 17, 18,       # Na, Mg, Al, (exclude Si=14), P, S, Cl, Ar
        19, 20,                              # K, Ca
        31,    32, 33, 34, 35, 36,            # Ga, (exclude Ge if needed? 32 is Ge and allowed here), As, Se, Br, Kr
        37, 38,                              # Rb, Sr
        49, 50, 51, 52, 53, 54,               # In, Sn, Sb, Te, I, Xe
        55, 56,                              # Cs, Ba
        81, 82, 83, 84, 86                   # Tl, Pb, Bi, Po, Rn
    }
    # Check that all atoms are allowed.
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atomic_numbers:
            return False, f"Atom with atomic number {atomic_num} is not allowed for a gas molecular entity"
    
    # Exclude molecules that contain silicon (14) or boron (5) even if they are in main–group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in {14, 5}:
            return False, f"Molecule contains element with atomic number {atom.GetAtomicNum()}, which is excluded"
    
    # Check molecular weight.
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight > 300:
        return False, f"Molecular weight too high for a gas ({mol_weight:.1f} Da > 300 Da)"
    
    # Check rotatable bonds.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds > 2:
        return False, f"Too many rotatable bonds ({rot_bonds}), indicating a flexible molecule unlikely to be a gas"
    
    # Check topological polar surface area (TPSA).
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    if mol_weight < 150:
        if tpsa >= 60:
            return False, f"TPSA too high ({tpsa:.1f} Å²) for a low-weight molecule"
    else:
        if tpsa >= 40:
            return False, f"TPSA too high ({tpsa:.1f} Å²), suggesting significant intermolecular interactions"
    
    # If the molecule contains any rings, only allow rings of size 3 (e.g. oxirane or cyclopropane).
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    if rings:
        for ring in rings:
            if len(ring) != 3:
                return False, "Molecule contains a ring larger than 3 atoms, which is unusual for a gas at STP"
    
    return True, "Molecule is a simple, low–weight, rigid main–group entity with low polarity, fitting expected characteristics of a gas at STP"

# Example usage (for testing purposes)
if __name__ == '__main__':
    test_smiles = [
       "[H][H]", "ClC=C", "CC(C)C", "FF", "C=C", "Cl[H]", "[H]\\C(C)=C(\\[H])C",
       "[Ne]", "[Kr]", "C(C(C(F)(F)F)(F)F)(C(F)(F)F)(F)F", "C1CO1",
       "[Xe]", "[O][O]", "[He]", "[Rn]", "[3H][3H]", "CC#C", "[220Rn]", "[C-]#[O+]",
       "[H]\\C(C)=C(/[H])C", "[Ar]", "O=C=O", "[222Rn]", "[6He]", "[C]",
       "FC=C", "CC", "CC(C)=C", "[1H][1H]", "[H]C([H])([H])[H]",
       "I[H]", "[O-][O+]=O", "ClCl"
    ]
    for smi in test_smiles:
        is_gas, reason = is_gas_molecular_entity(smi)
        print(f"SMILES: {smi:30s} -> {is_gas}, {reason}")
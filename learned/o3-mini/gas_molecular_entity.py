"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: gas molecular entity as defined by
"Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa)".

Heuristic changes from our previous approach:
  - Molecule must parse and be a single–fragment molecule.
  - We add explicit hydrogens so that hydride molecules (CH4, NH3, HCl, etc.) are not mistaken as single–atom entities.
  - Single–atom molecules (after adding H) are allowed only if they are noble gases.
  - All heavy atoms (non–hydrogen) must be from allowed main–group elements (we allow H explicitly), excluding Si and B.
  - Molecule must be small (molecular weight <=300 Da).
  - Molecule must be “rigid” (<=2 rotatable bonds).
  - The topological polar surface area (TPSA) must be low: if weight <150 Da then TPSA must be < 60 Å², else TPSA must be < 40 Å².
  - If the molecule contains any rings then only 3–membered rings are allowed.
  
This heuristic approach may still be imperfect.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string)
    qualifies as a gas molecular entity at STP.

    For our purposes, the molecule must:
      - Parse as a single–fragment molecule.
      - (If multi–atomic) contain only allowed main–group elements—with the proviso that
        compounds containing silicon (atomic number 14) or boron (atomic number 5) are excluded.
      - Have a molecular weight <=300 Da.
      - Have low topological polar surface area (TPSA): if MW <150, TPSA <60 Å²; else TPSA <40 Å².
      - Be “rigid” with <=2 rotatable bonds.
      - If rings are present, only rings of size 3 are allowed.
      
    Special note: Single–atom molecules are allowed only if that solitary atom is a noble gas.
    To avoid mis‐classifying hydride molecules (like CH4, NH3, HCl, HI), we add explicit hydrogens
    and then check if the molecule has only one atom.
    
    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
         bool: True if the molecule satisfies our criteria to be a gas molecular entity, False otherwise.
         str: Explanation for the decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure molecule is a single fragment (using the original mol)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (not a single molecular entity)"
        
    # Add explicit hydrogens to avoid mis‐interpreting hydride molecules as single–atom species.
    mol_with_H = Chem.AddHs(mol)
    
    # If, after adding hydrogens, the molecule still consists of a single atom,
    # then it is truly a single–atom molecular entity.
    if mol_with_H.GetNumAtoms() == 1:
        # Define noble gas atomic numbers: He (2), Ne (10), Ar (18), Kr (36), Xe (54), Rn (86)
        noble_gases = {2, 10, 18, 36, 54, 86}
        atom = mol_with_H.GetAtomWithIdx(0)
        if atom.GetAtomicNum() in noble_gases:
            return True, "Single noble gas atom: gas at STP"
        else:
            return False, f"Single atom of element with atomic number {atom.GetAtomicNum()} is not a gas at STP"
    
    # Check allowed heavy atoms: for atoms other than hydrogen, ensure they are in the allowed set.
    # Allowed main–group atomic numbers (including hydrogen) but we exclude Si (14) and B (5).
    allowed_atomic_numbers = {
         1,   2,    3, 4,   6, 7, 8, 9, 10,    # H, He, Li, Be, C, N, O, F, Ne
        11, 12,  13,    15, 16, 17, 18,        # Na, Mg, Al, P, S, Cl, Ar
        19, 20,                              # K, Ca
        31,    32, 33, 34, 35, 36,            # Ga, Ge, As, Se, Br, Kr
        37, 38,                              # Rb, Sr
        49, 50, 51, 52, 53, 54,               # In, Sn, Sb, Te, I, Xe
        55, 56,                              # Cs, Ba
        81, 82, 83, 84, 86                   # Tl, Pb, Bi, Po, Rn
    }
    # Use the molecule before adding H for checking elements (hydrogens are automatically allowed)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 1 and atomic_num not in allowed_atomic_numbers:
            return False, f"Atom with atomic number {atomic_num} is not allowed for a gas molecular entity"
        # Exclude silicon or boron if present (even if they are in main–group).
        if atomic_num in {14, 5}:
            return False, f"Molecule contains element with atomic number {atomic_num}, which is excluded"

    # Check molecular weight (using the original mol, where hydrogens are implicit)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight > 300:
        return False, f"Molecular weight too high for a gas ({mol_weight:.1f} Da > 300 Da)"
    
    # Check rotatable bonds (should be <=2 to ensure rigidity)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds > 2:
        return False, f"Too many rotatable bonds ({rot_bonds}), indicating a flexible molecule unlikely to be a gas"
    
    # Check the topological polar surface area (TPSA)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    if mol_weight < 150:
        if tpsa >= 60:
            return False, f"TPSA too high ({tpsa:.1f} Å²) for a low-weight molecule"
    else:
        if tpsa >= 40:
            return False, f"TPSA too high ({tpsa:.1f} Å²), suggesting significant intermolecular interactions"
    
    # Check rings: if any ring is present, only allow rings of size 3.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if rings:
        for ring in rings:
            if len(ring) != 3:
                return False, "Molecule contains a ring larger than 3 atoms, which is unusual for a gas at STP"
    
    return True, "Molecule is a simple, low–weight, rigid main–group entity with low polarity, fitting expected characteristics of a gas at STP"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
       "[H][H]", "ClC=C", "CC(C)C", "FF", "C=C", "Cl[H]", "[H]\\C(C)=C(\\[H])C",
       "[Ne]", "[Kr]", "C(C(C(F)(F)F)(F)F)(C(F)(F)F)(F)F", "C1CO1",
       "[Xe]", "[O][O]", "[He]", "[Rn]", "[3H][3H]", "CC#C", "[220Rn]", "[C-]#[O+]",
       "[H]\\C(C)=C(/[H])C", "[Ar]", "O=C=O", "[222Rn]", "[6He]", "[C]",
       "FC=C", "CC", "CC(C)=C", "[1H][1H]", "[H]C([H])([H])[H]",
       "I[H]", "[O-][O+]=O", "ClCl", "O=NN=O", "[H]OO[C]=O", "SCCC(=O)C",
       "C1(=CC1)C", "FC(F)(C(F)(F)F)C", "S(P(SC)(O)=O)C", "NC(=S)NCC=C",
       "C[As](C)([O-])=O", "CCOC(C)(C)C", "CC(C)[C@H](O)C(O)=O", "O=C=S",
       "CC(C(C)=O)C(O)=O", "OCCS", "ClC(Cl)=CC(O)=O", "O=C/C=C/C#CC#CC#CC",
       "OC(CCl)C(O)=O"
    ]
    
    for smi in test_smiles:
        is_gas, reason = is_gas_molecular_entity(smi)
        print(f"SMILES: {smi:40s} -> {is_gas:5s}, {reason}")
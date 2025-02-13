"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: gas molecular entity as defined by
"Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa)".

This implementation uses the following heuristic:
  1. The molecule’s SMILES must parse and represent a single fragment.
  2. All atoms must be from allowed main group elements.
  3. The molecule must have a low molecular weight (<=350 Da).
  4. The topological polar surface area (TPSA) must be low (<40 Å²) so that the molecule is unlikely to engage in extensive intermolecular hydrogen bonding.
  5. The molecule must have very few rotatable bonds (<=2), as gas molecules are generally rigid.
Note: This is a heuristic approach and might not capture all details.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string) qualifies as a gas molecular entity.
    
    A gas molecular entity (for our purposes) is assumed to be:
      - A single fragment.
      - Composed only of main group elements.
      - Low in molecular weight (<=350 Da).
      - Very low in polar surface area (<40 Å²).
      - Rigid (few rotatable bonds, <=2).
      
    Args:
        smiles (str): A SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria for a gas (gaseous at STP), False otherwise.
        str: Explanation for the classification decision.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule is a single connected entity.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (not a single molecular entity)"
    
    # Define allowed atomic numbers for main group elements.
    allowed_atomic_numbers = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10,    # H, He, Li, Be, B, C, N, O, F, Ne
        11, 12, 13, 14, 15, 16, 17, 18,     # Na, Mg, Al, Si, P, S, Cl, Ar
        19, 20,                           # K, Ca
        31, 32, 33, 34, 35, 36,             # Ga, Ge, As, Se, Br, Kr
        37, 38,                           # Rb, Sr
        49, 50, 51, 52, 53, 54,             # In, Sn, Sb, Te, I, Xe
        55, 56,                           # Cs, Ba
        81, 82, 83, 84, 86                # Tl, Pb, Bi, Po, Rn
    }
    
    # Check that all atoms in the molecule belong to the allowed main group elements.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Atom with atomic number {atom.GetAtomicNum()} is not from a main group element"

    # Calculate molecular weight.
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight > 350:
        return False, f"Molecular weight too high for a common gas ({mol_weight:.1f} Da > 350 Da)"
    
    # Calculate topological polar surface area (TPSA).
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    if tpsa >= 40:
        return False, f"TPSA too high ({tpsa:.1f} Å²), suggesting higher intermolecular interactions and lower volatility"
    
    # Calculate the number of rotatable bonds.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds > 2:
        return False, f"Too many rotatable bonds ({rot_bonds}), indicating a flexible molecule unlikely to be highly volatile"
    
    return True, "Molecule is a main group molecular entity with low weight, low polarity, and low flexibility, fitting expected characteristics of a gas at STP"
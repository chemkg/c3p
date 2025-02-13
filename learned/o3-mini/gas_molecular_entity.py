"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: gas molecular entity as defined by
"Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa)".

This implementation uses a simple heuristic:
  1. All atoms must be from main group elements (allowed atomic numbers defined below).
  2. The molecule must be a single fragment.
  3. The molecular weight must be below an arbitrary cutoff (350 Da) typical for low–molecular weight volatile compounds.
Note: This is a heuristic approach and might not capture all nuances.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule (given by its SMILES) is a gas molecular entity.

    Args:
        smiles (str): A SMILES string for the molecule.

    Returns:
        bool: True if the molecule meets our criteria for a gas molecular entity, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a single connected entity.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (not a single molecular entity)"

    # Define allowed atomic numbers for main group elements.
    # This set includes: H, He, Li, Be, B, C, N, O, F, Ne,
    # Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca,
    # and selected others common in main‐group molecules (e.g., Ga, Ge, As, Se, Br, Kr,
    # Rb, Sr, In, Sn, Sb, Te, I, Xe, Cs, Ba, Tl, Pb, Bi, Po, Rn).
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

    # Check each atom to ensure it belongs to a main group element.
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atomic_numbers:
            return False, f"Atom with atomic number {atomic_num} is not from a main group element"

    # Calculate the exact molecular weight using RDKit.
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    # Use an arbitrary threshold (350 Da) as a heuristic; most small gaseous molecules are below this value.
    if mol_weight > 350:
        return False, f"Molecular weight too high for common gases ({mol_weight:.1f} Da > 350 Da)"

    # If all conditions are met, we consider it as a gas molecular entity.
    return True, "Molecule is a main group molecular entity with low molecular weight and a single fragment, fitting expected characteristics of a gas at STP"
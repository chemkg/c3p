"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound (VOC)
Definition:
    Any organic compound having an initial boiling point less than or equal
    to 250 °C measured at a standard atmospheric pressure of 101.3 kPa.
    
Approach:
    Since predicting boiling point from the SMILES string is a challenging task
    and requires advanced quantitative models (e.g. using the Joback method),
    here we use a simple heuristic based on molecular weight:
      - We first ensure that the molecule is organic (i.e. it contains at least one carbon atom).
      - Then, we estimate its molecular weight.
      - We assume that if the molecular weight is <= 250 Da the compound is likely volatile.
    Note: This is only a rough screen and may not be reliable for every compound.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    
    A VOC is defined as any organic compound having an initial boiling point <= 250°C 
    measured at 101.3 kPa. Here we use molecular weight as a crude proxy.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is predicted to be a volatile organic compound, False otherwise.
        str: Reason for classification
    """
    
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check that the molecule is organic (contains at least one carbon atom)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found; not considered an organic compound"
    
    # Calculate the molecular weight (exact, in Daltons)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Use a heuristic: many small organic compounds (with lower molecular weight)
    # are volatile. This cutoff is crude; many factors affect boiling point.
    threshold = 250.0  # Dalton cutoff chosen based on the approximate definition.
    if mol_wt <= threshold:
        return True, f"Molecular weight ({mol_wt:.1f} Da) is below threshold; likely has an initial boiling point <= 250°C"
    else:
        return False, f"Molecular weight ({mol_wt:.1f} Da) exceeds threshold; likely has an initial boiling point > 250°C"
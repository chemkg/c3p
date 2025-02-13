"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule 
Definition: “A macromolecule is a molecule of high relative molecular mass, the structure of which essentially comprises the multiple repetition of units derived, actually or conceptually, from molecules of low relative molecular mass.”
Heuristic:
  • The SMILES is parsed successfully.
  • Its exact molecular weight is above 1000 Da.
  • It contains a high number (≥ 50) of heavy (non‐hydrogen) atoms.
If these conditions hold, we consider the molecule to be a macromolecule.
Note: This is only a heuristic approach.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) is a macromolecule.
    Our heuristic requires that:
      - The molecule is valid.
      - It has a high relative molecular weight (≥ 1000 Da).
      - It has a relatively large number of heavy atoms (≥ 50).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule passes the heuristic for a macromolecule, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the input SMILES string to a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Calculate the exact molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Count heavy (non-hydrogen) atoms.
    heavy_atoms = mol.GetNumHeavyAtoms()
    
    # Check if the molecular weight is high enough.
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da); macromolecules are usually ≥ 1000 Da."
    
    # Check if the molecule contains enough heavy atoms.
    if heavy_atoms < 50:
        return False, f"Not enough heavy atoms ({heavy_atoms}); macromolecules typically contain many atoms."
    
    # Optionally, one might try to detect repeating substructural patterns.
    # However, such an analysis is complex and beyond the scope of this simple heuristic.
    
    return True, f"Molecule qualifies as a macromolecule (Molecular weight: {mol_wt:.2f} Da; Heavy atoms: {heavy_atoms})."
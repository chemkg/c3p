"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin contains a skeleton of four pyrrole nuclei united through the alpha-positions 
    by four methine (-CH=) groups that form a macrocyclic structure, often coordinating a metal ion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved Porphyrin SMARTS pattern
    # Pattern description: Four pyrrole rings (c1[nH]cc1) with methine bridges (-C=) connecting them,
    # typically forming a central cavity which might contain a metal ion.
    porphyrin_pattern = Chem.MolFromSmarts("n1c(ccc1)-c2c(cnc3[nH]c(cc23)C)-c4[nH]ccc5cnc(cc5-c6ncccc6)c4")
    
    # Check for porphyrin macrocycle structure
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "Structure does not match porphyrin macrocycle"

    # Check for potential metal coordination, allowing for common metals found in porphyrins
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in ['Fe', 'Mg', 'Zn', 'Co', 'Ni', 'Cu', 'Mn']]
    if len(metal_atoms) == 0:
        return False, "No metal ion commonly associated with porphyrins found"

    return True, "Contains porphyrin macrocycle with potential metal coordination"
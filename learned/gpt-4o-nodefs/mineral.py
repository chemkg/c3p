"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Minerals typically have metal cations paired with inorganic anions,
    often with complex ion structures and/or hydration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expand list of metal elements typically found in minerals
    metals = set(['Ca', 'Mg', 'Fe', 'Ni', 'Ba', 'Zn', 'Cu', 'K', 
                  'Cs', 'Na', 'Al', 'Sb', 'La', 'Pb', 'Si', 'As', 'Pd', 'V'])
    metal_atoms_in_mol = {atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() in metals}
    
    if not metal_atoms_in_mol:
        return False, "No typical metal elements found"
    
    # Look for common inorganic anions or elements like S, Cl, F, and OH groups
    patterns = [
        '[O-]S([O-])(=O)=O',  # sulfate
        '[O-]P([O-])([O-])=O',  # phosphate
        '[O-]C([O-])=O',      # carbonate
        'Cl',                 # chloride
        'F',                  # fluoride
        '[S-]',               # sulfide used in pyrite etc.
        '[O-]',               # oxide
        'N',                  # for nitrides, cyanides etc.
        '[B-]',               # borates
        '~[Si]',              # silicates
    ]
    anion_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in patterns)
    
    if not anion_found:
        return False, "No typical inorganic anions or complex elements found"
    
    # Consider the presence of non-hydrocarbon organic structure as a warning sign
    organic_elements = set(['C', 'H'])
    organic_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in organic_elements]
    
    if len(organic_atoms) > 5:  # Arbitrarily assuming most minerals have minimal organic components
        organic_structure = Chem.MolFromSmiles('C')  # Simplistic pattern to detect organic content
        if any(mol.HasSubstructMatch(organic_structure) for organic_structure in organic_atoms):
            return False, "Too many organic components indicative of a non-mineral"

    # Check for representative hydration patterns (water molecules)
    hydrate_patterns = [
        Chem.MolFromSmarts("O"),
    ]
    hydrate_found = any(mol.HasSubstructMatch(hydrate_pattern) for hydrate_pattern in hydrate_patterns)

    return True, f"Contains metal elements ({', '.join(sorted(metal_atoms_in_mol))}) and inorganic anions{' with hydration' if hydrate_found else ''}"
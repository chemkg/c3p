"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:50828 aromatic primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is a primary alcohol in which the alcoholic hydroxy group
    is attached to a carbon which is itself bonded to an aromatic ring or aromatic heterocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and remove tautomers
    mol = Chem.MolFromSmiles(smiles, removeHs=False)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.RemoveHs(mol)
    
    # Look for primary alcohol (-CH2OH) group
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2;X3][OX2H]")
    primary_alcohol_matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    if not primary_alcohol_matches:
        return False, "No primary alcohol group found"
    
    # Check if any primary alcohol carbon is attached to an aromatic ring or heterocycle
    for match in primary_alcohol_matches:
        atom_idx = match[0]  # Index of primary alcohol carbon
        atom = mol.GetAtomWithIdx(atom_idx)
        if any(mol.GetAtomWithIdx(nbr_idx).GetIsAromatic() for nbr_idx in atom.GetNeighbors()):
            return True, "Contains a primary alcohol attached to an aromatic ring or heterocycle"
        
        # Check for aromatic heterocycles
        for nbr_idx in atom.GetNeighbors():
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            if nbr_atom.GetIsAromatic() and nbr_atom.GetAtomicNum() in [5, 6, 7, 8]:  # N, C, O, S
                return True, "Contains a primary alcohol attached to an aromatic heterocycle"
    
    return False, "Primary alcohol not attached to an aromatic ring or heterocycle"

# Example usage:
print(is_aromatic_primary_alcohol("Nc1nc(=O)[nH]cc1CO"))  # True, "Contains a primary alcohol attached to an aromatic heterocycle"
print(is_aromatic_primary_alcohol("CCCCCC"))  # False, "No primary alcohol group found"
print(is_aromatic_primary_alcohol("OC(=O)c1ccccc1CO"))  # False, "Primary alcohol not attached to an aromatic ring or heterocycle"
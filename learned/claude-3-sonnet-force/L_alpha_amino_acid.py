"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:32567 L-alpha-amino acid
Any alpha-amino acid having L-configuration at the alpha-carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alpha-amino acid pattern (C-N-C(=O)-O)
    amino_acid_pattern = Chem.MolFromSmarts("[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No alpha-amino acid pattern found"

    # Find chiral centers in the molecule
    chiral_centers = Chem.FindMolChiralCenters(mol)

    # Identify the alpha-carbon
    alpha_carbon = None
    for center in chiral_centers:
        atom = mol.GetAtomWithIdx(center)
        if atom.GetSymbol() == "C" and len(atom.GetNeighbors()) == 3:
            alpha_carbon = atom
            break

    if alpha_carbon is None:
        return False, "No alpha-carbon found"

    # Determine chirality at alpha-carbon
    try:
        conf = alpha_carbon.GetProp("_CIPCode")
        if conf == "R":
            return True, "L-configuration at alpha-carbon"
        elif conf == "S":
            return False, "D-configuration at alpha-carbon"
        else:
            return False, "Unable to determine chirality at alpha-carbon"
    except KeyError:
        return False, "Unable to determine chirality at alpha-carbon"
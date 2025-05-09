"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid has both an amino group (-NH2) and a carboxyl group (-COOH)
    but is not one of the standard 20 proteinogenic amino acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of both amino and carboxyl groups
    amino_pattern = Chem.MolFromSmarts("[NH2,NH1;!$(NC=O)]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # More comprehensive list of proteinogenic amino acids
    proteinogenic_amino_acids = [
        # Standard 20 amino acids
        "C[C@@H](C(=O)O)N",  # Alanine
        "N[C@@H](CC(=O)O)C(=O)O",  # Aspartic acid
        "N[C@@H](CCCNC(=N)N)C(=O)O",  # Arginine
        "N[C@@H](CC(=O)N)C(=O)O",  # Asparagine
        "N[C@@H](CS)C(=O)O",  # Cysteine
        "N[C@@H](CCC(=O)O)C(=O)O",  # Glutamic acid
        "N[C@@H](CCCNC(=O)N)C(=O)O",  # Glutamine
        "N[C@@H](CC1=CNC=N1)C(=O)O",  # Histidine
        "N[C@@H](CC(C)C)C(=O)O",  # Isoleucine
        "N[C@@H](CCCCN)C(=O)O",  # Lysine
        "N[C@@H](CCSC)C(=O)O",  # Methionine
        "N[C@@H](Cc1ccccc1)C(=O)O",  # Phenylalanine
        "N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O",  # Tryptophan
        "N[C@@H](CC1=CC=C(O)C=C1)C(=O)O",  # Tyrosine
        "N[C@@H](C)C(=O)O",  # Glycine
        "N[C@@H](CC(C)(C)C)C(=O)O",  # Valine
        "N[C@@H](CC(O)C)C(=O)O",  # Threonine
        "N[C@@H](CC1=CC=C(O)C=C1)C(=O)O",  # Tyrosine (alternative)
        "N[C@@H](CC1=CNC=N1)C(=O)O",  # Histidine (alternative)
        "N[C@@H](CO)C(=O)O",  # Serine
        "N[C@@H](CC(=O)O)C(=O)O",  # Aspartic acid (alternative)
        "N[C@@H](CCCNC(=O)N)C(=O)O",  # Glutamine (alternative)
        "N[C@@H](CC1=CC=CC=C1)C(=O)O",  # Phenylalanine (alternative)
        "N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O",  # Tryptophan (alternative)
    ]

    # Check against proteinogenic patterns using more precise matching
    for paa_smiles in proteinogenic_amino_acids:
        paa_mol = Chem.MolFromSmiles(paa_smiles)
        if mol.HasSubstructMatch(paa_mol):
            # Additional check: if the match is exact, it's proteinogenic
            if Chem.MolToSmiles(mol) == Chem.MolToSmiles(paa_mol):
                return False, "Molecule is a standard proteinogenic amino acid"
            # If it's a substructure match but not exact, it might be a modified version
            # so we continue checking

    # More specific patterns for non-proteinogenic features
    non_standard_patterns = [
        "[OH]",  # Hydroxyl group
        "[P](=O)(O)(O)",  # Phosphate group
        "[S](=O)(=O)",  # Sulfoxide/sulfone
        "[N+](=O)[O-]",  # Nitro group
        "[C]=[C]",  # Double bonds in side chains
        "[C]#[C]",  # Triple bonds in side chains
        "[C]=O",  # Ketones
        "[C](=O)O",  # Carboxylates
        "[C](=O)N",  # Amides
        "[C](=O)S",  # Thioesters
        "[N](=O)",  # Nitroso
        "[C](=O)[C]",  # Ketones in side chains
        "[C](=O)[N]",  # Amides in side chains
        "[C](=O)[O]",  # Esters
        "[C](=O)[S]",  # Thioesters
        "[C](=O)[C](=O)",  # Dicarbonyls
        "[C](=O)[C](=O)[C]",  # Tricarbonyls
        "[C](=O)[C](=O)[C](=O)",  # Tetracarbonyls
        "[C](=O)[C](=O)[C](=O)[C]",  # Pentacarbonyls
        "[C](=O)[C](=O)[C](=O)[C](=O)",  # Hexacarbonyls
    ]

    # Check for non-standard features
    for pattern in non_standard_patterns:
        patt_mol = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt_mol):
            return True, "Contains non-standard functional groups"

    # If no non-standard features are found, it might still be a non-proteinogenic amino acid
    # but with a less obvious modification
    return True, "Contains amino and carboxyl groups but is not a standard proteinogenic amino acid"
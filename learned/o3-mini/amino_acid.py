"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: Amino acids
Definition: A carboxylic acid containing one or more amino groups.
Examples include standard α–amino acids as well as non‐canonical structures.
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a molecule that has a carboxylic acid group
    and one or more amino (free amine) groups (i.e. not part of an amide).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an amino acid, False otherwise.
        str: A reason for the classification.
    """
    # Parse SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for carboxylic acid.
    # The first pattern matches the protonated carboxyl group: C(=O)OH.
    # The second pattern matches the deprotonated carboxylate: C(=O)[O-].
    acid_smarts1 = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_smarts2 = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    # Check if the molecule contains a carboxylic acid group.
    has_acid = mol.HasSubstructMatch(acid_smarts1) or mol.HasSubstructMatch(acid_smarts2)
    if not has_acid:
        return False, "No carboxylic acid group found"
    
    # Define a SMARTS pattern for an amino group.
    # The general amine group (sp3-hybridized nitrogen) is represented by [NX3].
    # We exclude cases where the nitrogen is directly bonded to a carbonyl (N-C(=O))
    # i.e. where it forms an amide, using !$(N[C](=O)).
    amine_smarts = Chem.MolFromSmarts("[NX3;!$(N[C](=O))]")
    
    # Check if at least one amino group is found.
    if not mol.HasSubstructMatch(amine_smarts):
        return False, "No free amino group found"
    
    return True, "Molecule contains a carboxylic acid and at least one free amino group"
"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl Sulfide
Definition: Any aliphatic sulfide in which at least one of the organyl groups attached 
to the sulfur is a methyl group.
Examples include S-methyl-L-ergothioneine, 1-(methylthio)ribulose 5-phosphate, and others.
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a given molecule is a methyl sulfide.
    A methyl sulfide is defined as an aliphatic (non‐aromatic) sulfide in which the sulfur 
    atom is bonded to two carbon groups and at least one of these is a methyl group (i.e. 
    a carbon atom with three hydrogens attached).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a methyl sulfide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are present in the structure.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern that matches an aliphatic sulfur (non‐aromatic, 2-connected)
    # bonded to a methyl group (a tetrahedral carbon with exactly three hydrogens)
    # and any other carbon.
    # The pattern [C;X4&H3] matches a carbon that is sp3 and has exactly 3 attached hydrogens,
    # [S;!a;X2] matches a non-aromatic sulfur with exactly two bonds,
    # and [#6] matches any carbon atom.
    pattern = Chem.MolFromSmarts("[C;X4&H3]-[S;!a;X2]-[#6]")
    
    # Check if the pattern is found in the molecule.
    if mol.HasSubstructMatch(pattern):
        return True, "Found an aliphatic sulfide with at least one methyl group attached to sulfur."
    
    return False, "No suitable aliphatic sulfide with a methyl group attached was found."
"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A fatty acid is here defined as a molecule containing a carboxylic acid group and a largely aliphatic (carbon-rich) region.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid, defined as a fatty acid containing at least one ring anywhere in its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a cyclic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group.
    # This simple SMARTS pattern looks for a carbonyl directly attached to a hydroxyl.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group detected; not a fatty acid"
    
    # Count carbon atoms to ensure the molecule has a sufficient "fatty" (long-chain) character.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Setting an arbitrary threshold of 7 carbons to differentiate typical fatty acids from smaller acids.
    if c_count < 7:
        return False, f"Too few carbons (found {c_count}) for a fatty acid"
    
    # Check if the molecule contains at least one ring.
    if mol.GetRingInfo().NumRings() < 1:
        return False, "No ring detected in the structure; not a cyclic fatty acid"
    
    # If all checks pass, classify the molecule as a cyclic fatty acid.
    return True, "Contains carboxylic acid group and ring structure; qualifies as a cyclic fatty acid"
"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition: Any fatty acid containing at least one C=C or C#C bond.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is defined here as a fatty acid (a molecule containing a carboxyl group)
    that has at least one carbon-carbon double or triple bond.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a carboxylic acid group:
    # We try two patterns that commonly represent acid groups: neutral and deprotonated.
    acid_pattern_neutral = Chem.MolFromSmarts("C(=O)O")
    acid_pattern_deprotonated = Chem.MolFromSmarts("C(=O)[O-]")
    
    has_acid = mol.HasSubstructMatch(acid_pattern_neutral) or mol.HasSubstructMatch(acid_pattern_deprotonated)
    if not has_acid:
        return False, "Molecule does not contain a carboxylic acid group, so not a fatty acid"
    
    # Check for unsaturation: we need at least one carbon-carbon double or triple bond.
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    triple_bond_pattern = Chem.MolFromSmarts("[#6]#[#6]")
    
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    has_triple_bond = mol.HasSubstructMatch(triple_bond_pattern)
    
    if has_double_bond or has_triple_bond:
        # Provide details on which unsaturation was found.
        if has_double_bond and has_triple_bond:
            unsat_detail = "carbon-carbon double and triple bond(s)"
        elif has_double_bond:
            unsat_detail = "carbon-carbon double bond(s)"
        else:
            unsat_detail = "carbon-carbon triple bond(s)"
        return True, f"Molecule is a fatty acid and contains unsaturation: {unsat_detail}"
    else:
        return False, "Fatty acid, but does not contain any carbon-carbon double or triple bonds (unsaturation)"
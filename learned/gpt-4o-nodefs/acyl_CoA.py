"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str) -> (bool, str):
    """
    Classifies the molecule as acyl-CoA based on SMILES string.
    Acyl-CoA is characterized by a fatty acid chain bound to Coenzyme A via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as acyl-CoA, False otherwise.
        str: Reason for the classification.
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated SMARTS pattern to identify key acyl-CoA components
    # Focus on the CoA structure: adenine, ribose, phosphate, and pantetheine parts
    coa_pattern = Chem.MolFromSmarts(
        "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](COP(=O)(O)O)[C@H](O)[C@@H]1OP(=O)(O)O[C@@H]1N=C(N)N=C2C=CN(C2=N1)C"
    )
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not identified"

    # Update thioester pattern to confirm acyl linkage with flexibility in structures
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester bond with CoA moiety not identified"

    return True, "Contains CoA moiety and thioester bond, characteristic of acyl-CoA"
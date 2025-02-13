"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies chemical entities of the class carbamate ester.
Definition: any ester of carbamic acid or its N-substituted derivatives.
Characteristic substructure: R–O–C(=O)–NR′.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as any ester of carbamic acid or its N-substituted derivatives,
    containing the substructure R-O-C(=O)-NR' (the O-C(=O)-N motif).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a carbamate ester functional group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carbamate ester: O-C(=O)-N
    # This pattern looks for an oxygen atom directly bonded to a carbonyl carbon,
    # which in turn is bonded to a nitrogen atom.
    carbamate_smarts = "O-C(=O)N"
    carbamate_pattern = Chem.MolFromSmarts(carbamate_smarts)
    if carbamate_pattern is None:
        return False, "Error creating SMARTS pattern for carbamate ester"

    # Check if the molecule contains the carbamate ester substructure.
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Molecule contains a carbamate ester group (O-C(=O)-N)."
    else:
        return False, "Carbamate ester group (O-C(=O)-N) not found in the molecule."
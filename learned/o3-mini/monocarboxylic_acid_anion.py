"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: Monocarboxylic acid anion 
Definition: A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
A genuine monocarboxylic acid anion should contain exactly one deprotonated carboxylate group ([O-]C(=O) attached to a carbon)
and the overall molecule should have a net charge of -1. Also, the molecule should be a single fragment (i.e. not a salt).
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    Requirements:
      - The molecule should contain exactly one deprotonated carboxylate group.
      - The overall net formal charge should be -1.
      - The molecule should be a single fragment (i.e. not include counterions).
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a monocarboxylic acid anion, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule is a single fragment.
    # If there are multiple fragments (indicated by a dot in the SMILES),
    # then likely the molecule is a salt or a mixture.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule contains multiple fragments; likely a salt rather than a bare anion"

    # Compute the net formal charge of the molecule.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -1:
        return False, f"Molecule net charge is {net_charge}; expected -1 for a monocarboxylic acid anion"

    # Define the SMARTS pattern for a deprotonated carboxylate:
    # [O-] attached to a carbon which has a double bond to an oxygen (i.e. the C=O group)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)[#6]")
    if carboxylate_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches of the carboxylate group.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    num_matches = len(matches)
    
    # There must be exactly one deprotonated carboxylate group.
    if num_matches == 0:
        return False, "No deprotonated carboxylate (COO-) group found"
    elif num_matches > 1:
        return False, f"Found {num_matches} carboxylate groups; expected exactly one for a monocarboxylic acid anion"
    
    # Passed all the checks: exactly one carboxylate group in a single-fragment molecule with net charge -1.
    return True, "Contains exactly one deprotonated carboxylate (COO-) group and net charge -1, consistent with a monocarboxylic acid anion"

# Example usage (can be removed or commented out when integrating into a larger system):
if __name__ == "__main__":
    # Test with a known monocarboxylic acid anion SMILES.
    test_smiles = "CCCC(O)CC([O-])=O"  # 3-hydroxyhexanoate, expected True
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(f"Result: {result}\nReason: {reason}")
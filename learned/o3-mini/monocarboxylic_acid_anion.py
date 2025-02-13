"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
A genuine monocarboxylic acid anion should contain exactly one deprotonated carboxylate group ([O-]C(=O)–R)
and the overall molecule should have a net formal charge of -1. In addition, the molecule should be a single fragment,
and it should not contain any positively charged atoms that would hint at internal salt formation.
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    The requirements are:
      - The molecule must be parsed successfully.
      - The molecule should be a single fragment (i.e. not include counterions).
      - None of the atoms should carry a positive formal charge.
      - The molecule should have a net formal charge of -1.
      - There must be exactly one deprotonated carboxylate group present.
        The carboxylate group is defined as an oxygen with a -1 charge directly bonded to a carbon,
        where that carbon is double-bonded to another oxygen.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule qualifies as a monocarboxylic acid anion, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule contains multiple fragments; likely a salt rather than a bare anion"
    
    # Check that no atom has a positive formal charge.
    if any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "Molecule contains positively charged atom(s), not consistent with a bare monocarboxylic acid anion"
    
    # Compute the net formal charge of the molecule.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -1:
        return False, f"Molecule net charge is {net_charge}; expected -1 for a monocarboxylic acid anion"
    
    # Define the SMARTS pattern for a deprotonated carboxylate.
    # This pattern matches an oxygen with -1 charge directly bonded to a carbon that has a double bond to a second oxygen.
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)[#6]")
    if carboxylate_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches of the carboxylate group.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    num_matches = len(matches)
    if num_matches == 0:
        return False, "No deprotonated carboxylate (COO-) group found"
    elif num_matches > 1:
        return False, f"Found {num_matches} carboxylate groups; expected exactly one for a monocarboxylic acid anion"
    
    # Verify that the match truly represents a carboxylate functionality.
    # Our SMARTS "[O-]C(=O)[#6]" is expected to return a tuple of indices:
    #   index0: oxygen with -1 charge,
    #   index1: the carboxyl carbon,
    #   index2: the carbon substituent attached to that carboxyl carbon.
    match = matches[0]
    o_minus_idx, c_idx, r_idx = match
    o_minus_atom = mol.GetAtomWithIdx(o_minus_idx)
    c_atom = mol.GetAtomWithIdx(c_idx)
    
    # Confirm that the oxygen is indeed deprotonated.
    if o_minus_atom.GetAtomicNum() != 8 or o_minus_atom.GetFormalCharge() != -1:
        return False, "The matched oxygen is not a deprotonated oxygen (O-)"
    
    # Check that the carboxyl carbon has a double-bonded oxygen neighbor.
    found_double_bond = False
    for nbr in c_atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
        if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 1:
            found_double_bond = True
            break
    if not found_double_bond:
        return False, "Carboxylate carbon is not attached via a double bond to an oxygen (missing C=O)"
    
    # All criteria passed.
    return True, "Contains exactly one deprotonated carboxylate (COO-) group, no interfering positive charges, and net charge -1 consistent with a monocarboxylic acid anion"

# Example usage (this section can be removed or commented out when integrating into a larger system):
if __name__ == "__main__":
    # Test with a known monocarboxylic acid anion SMILES (e.g., 3-hydroxyhexanoate).
    test_smiles = "CCCC(O)CC([O-])=O"
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(f"Test SMILES: {test_smiles}\nResult: {result}\nReason: {reason}")
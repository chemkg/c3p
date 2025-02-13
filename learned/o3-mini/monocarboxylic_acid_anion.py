"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxy group of a monocarboxylic acid is deprotonated.
A genuine monocarboxylic acid anion should contain exactly one deprotonated carboxylate group (R-C(=O)[O-])
and have a net formal charge of -1. The molecule should be a single fragment and must not contain any positively charged atoms.
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    
    The requirements are:
      - The SMILES string is parseable.
      - The molecule is a single fragment (i.e. no separate counterions).
      - None of the atoms have a positive formal charge.
      - The net formal charge of the molecule is -1.
      - There is exactly one deprotonated carboxylate group (R-C(=O)[O-]) present.
        We define this functional group using an explicit SMARTS with mapping:
          [O-:1]-[C:2](=[O:3])
        where:
          :1 is the deprotonated oxygen (O with -1 charge),
          :2 is the carboxyl carbon,
          :3 is the carbonyl oxygen (C=O).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a monocarboxylic acid anion, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is a single fragment (no counterions or salts)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule contains multiple fragments; likely a salt rather than a bare anion"
    
    # Check that no atom has a positive formal charge.
    if any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "Molecule contains one or more positively charged atoms"
    
    # Calculate the net formal charge; should be exactly -1.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -1:
        return False, f"Molecule net charge is {net_charge}; expected -1 for a monocarboxylic acid anion"
    
    # Create a SMARTS pattern with explicit atom mapping to capture exactly three atoms:
    #   [O-:1] is the deprotonated oxygen (should have a -1 formal charge),
    #   [C:2] is the carboxyl carbon,
    #   [O:3] is the carbonyl oxygen connected by a double bond.
    carboxylate_pattern = Chem.MolFromSmarts("[O-:1]-[C:2](=[O:3])")
    if carboxylate_pattern is None:
        return False, "Error creating SMARTS pattern for the carboxylate group"
    
    # Find all substructure matches of the carboxylate group.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    num_matches = len(matches)
    if num_matches == 0:
        return False, "No deprotonated carboxylate group found"
    elif num_matches > 1:
        return False, f"Found {num_matches} carboxylate groups; expected exactly one"
    
    # Now unpack the only match; with our explicit mapping it should exactly contain three indices.
    match = matches[0]
    if len(match) != 3:
        return False, f"Unexpected number of atoms in the carboxylate match: {len(match)} (expected 3)"
    o_minus_idx, c_idx, dobly_o_idx = match
    o_minus_atom = mol.GetAtomWithIdx(o_minus_idx)
    c_atom = mol.GetAtomWithIdx(c_idx)
    dobly_o_atom = mol.GetAtomWithIdx(dobly_o_idx)
    
    # Confirm that the deprotonated oxygen has the expected -1 charge.
    if o_minus_atom.GetAtomicNum() != 8 or o_minus_atom.GetFormalCharge() != -1:
        return False, "The matched oxygen does not have a -1 formal charge as required"
    
    # Verify that the carboxyl carbon has a double-bond to the second oxygen.
    found_double_bond = False
    for bond in c_atom.GetBonds():
        # Check if bond is double and the bonded neighbor is the carbonyl oxygen.
        if bond.GetBondTypeAsDouble() == 2.0:
            neighbor = bond.GetOtherAtom(c_atom)
            if neighbor.GetIdx() == dobly_o_atom.GetIdx():
                found_double_bond = True
                break
    if not found_double_bond:
        return False, "Carboxyl carbon is not double-bonded to the expected oxygen (C=O)"
    
    return True, "Contains exactly one deprotonated carboxylate group with net charge -1 and meets all criteria"

# Example usage:
if __name__ == "__main__":
    # Test with a sample monocarboxylic acid anion, e.g. 3-hydroxyhexanoate.
    test_smiles = "CCCC(O)CC([O-])=O"
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(f"Test SMILES: {test_smiles}\nResult: {result}\nReason: {reason}")
"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo group is located at the 2-position.
That is, the carboxylate group (C(=O)[O-]) is directly attached to an alpha carbon which bears a carbonyl (=O) substituent.
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    It searches for a carboxylate group (C(=O)[O-]) that is attached to an alpha-carbon bearing an additional carbonyl,
    which corresponds to the oxo group at the 2-position.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule fits the definition, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carboxylate anion group: C(=O)[O-]
    carboxylate_smarts = "C(=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    
    # Find all matches for the carboxylate group
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not matches:
        return False, "No carboxylate anion (C(=O)[O-]) group found"

    # For each carboxylate group found, attempt to find an alpha carbon that bears an additional carbonyl (=O) substituent.
    for match in matches:
        # match[0] corresponds to the carboxylate carbon in the SMARTS "C(=O)[O-]"
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Look at neighbors of the acid carbon to get the alpha carbon (exclude the oxygens that are part of the carboxylate group)
        for neighbor in acid_carbon.GetNeighbors():
            # The alpha carbon usually is carbon (atomic number 6)
            if neighbor.GetAtomicNum() != 6:
                continue
            alpha_carbon = neighbor
            # Check if this alpha carbon has a double bond to an oxygen (the oxo group).
            # We iterate over its bonds.
            has_oxo = False
            for bond in alpha_carbon.GetBonds():
                # We require a double bond:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other_atom = bond.GetOtherAtom(alpha_carbon)
                    # Check that the double-bonded partner is an oxygen
                    if other_atom.GetAtomicNum() == 8:
                        # We need to make sure that this oxygen is not the one already in the acid group.
                        # Since the acid_carbon already is double-bonded to one oxygen in the carboxylate, we want a separate one.
                        # We can compare the indices.
                        if other_atom.GetIdx() not in match[1:]:
                            has_oxo = True
                            break
            if has_oxo:
                return True, ("Found carboxylate group attached to an alpha carbon that bears a "
                              "double-bonded oxygen (oxo) substituent, fitting the 2-oxo monocarboxylic acid anion definition.")
    
    return False, "No alpha carbon with an additional oxo substituent (at the 2-position) found adjacent to a carboxylate group"

# Example usage (you can test with one of the provided SMILES):
if __name__ == "__main__":
    test_smiles = "CC(=O)CC(=O)C([O-])=O"  # acetylpyruvate as an example
    result, reason = is_2_oxo_monocarboxylic_acid_anion(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)
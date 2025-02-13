"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol
Definition: Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.
They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.
"""

from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal in which an intramolecular addition of a hydroxy group
    to a carbonyl leads to a 1-oxacycloalkan-2-ol motif (or an unsaturated analogue).

    The detection uses a SMARTS pattern that looks for a ring carbon which is not a carbonyl,
    having an exocyclic hydroxyl group bonded (i.e. [OX2H] that is not in a ring) and also bonded 
    to an oxygen atom that is part of the ring. This pattern is a simple heuristic and may not cover 
    all cases in complex molecules.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a lactol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for a lactol center:
    # - [C;R;!$(C=O)] describes a ring carbon that is not a carbonyl.
    # - ([OX2H;!R]) is an exocyclic hydroxyl group (oxygen with H not in a ring).
    # - [O;R] is a ring oxygen, making the -OR portion of the cyclic hemiacetal.
    lactol_smarts = "[C;R;!$(C=O)]([OX2H;!R])[O;R]"
    lactol_pat = Chem.MolFromSmarts(lactol_smarts)
    if lactol_pat is None:
        return False, "Error in lactol SMARTS pattern"

    # Check if the molecule has at least one lactol-like substructure match
    if mol.HasSubstructMatch(lactol_pat):
        return True, "Found lactol moiety pattern (cyclic hemiacetal) in the molecule"
    else:
        return False, "No lactol moiety pattern (cyclic hemiacetal) detected in the molecule"

# Example usage:
if __name__ == "__main__":
    # Test a few known lactol SMILES (e.g., a sugar in its cyclic hemiacetal form)
    test_smiles = [
        "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-D-fructopyranose (a cyclic hemiacetal)
        "O=C1OC(CO)C(O)C1"  # an example that may or may not be a lactol based on ring closure
    ]
    for s in test_smiles:
        result, reason = is_lactol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")
"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines
Definition: Any piperazinone that has a piperazine-2,5-dione skeleton.
This implementation attempts to detect a six-membered ring having exactly
two nitrogen atoms and two carbonyl groups in positions corresponding to a cyclic dipeptide.
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines whether a given molecule (via SMILES) contains a piperazine-2,5-dione skeleton.
    The improved approach uses a refined SMARTS pattern to look for a six-membered ring (cyclic dipeptide core)
    with two nitrogen atoms and two C=O groups in the appropriate positions. RDKit is directed
    to ignore differences in stereochemistry for a more inclusive match.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if a piperazine-2,5-dione skeleton is detected, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define an improved SMARTS pattern for a 2,5-diketopiperazine core.
    # Explanation:
    #  - The pattern uses ring closure (1) to insist on a single six-membered ring.
    #  - [NX3;R] designates a (trivalent) nitrogen in a ring.
    #  - [CX3](=O) designates an sp2-hybridized carbon that is double-bonded to oxygen (a carbonyl).
    #  - [CX4;R] designates an sp3-hybridized carbon in the ring.
    # The intended order (cycling around the ring) is:
    #    N - C(=O) - C - N - C(=O) - C  (and then back to the first N).
    diketopiperazine_smarts = "[NX3;R]1[CX3](=O)[CX4;R][NX3;R][CX4;R][CX3](=O)[CX4;R]1"
    
    pattern = Chem.MolFromSmarts(diketopiperazine_smarts)
    if pattern is None:
        return False, "Error in constructing the SMARTS."

    # Use a substructure search ignoring stereochemistry.
    matches = mol.GetSubstructMatches(pattern, useChirality=False)
    if matches:
        # It is possible that a false positive match may occur on a fragment that is not the genuine core.
        # (For example, extra substituents might yield a match even if the diketopiperazine moiety is not the main scaffold.)
        # Here we check that at least one match involves six atoms (the six atoms of the core ring).
        for match in matches:
            if len(match) == 6:
                return True, "Molecule contains a piperazine-2,5-dione skeleton."
        # If no match with 6 atoms was found, then our pattern may have only matched a partial fragment.
        return False, "Substructure found, but it does not form a complete 6-membered diketopiperazine ring."
    else:
        return False, "Molecule does not contain a piperazine-2,5-dione skeleton."

# Example usage:
if __name__ == "__main__":
    # You can test with the provided SMILES strings.
    test_smiles = [
        ("Brocazine F", "S1S[C@]23N([C@@H]4[C@@H](O)C=C[C@@H]([C@H]4C2)O)C([C@]15N([C@@H]6[C@@H](O)C=CC([C@H]6C5)=O)C3=O)=O"),
        ("piperazine-2,5-dione", "O=C1CNC(=O)CN1"),
        ("(S)-Mephenytoin", "O=C1N(C(=O)N[C@@]1(CC)C2=CC=CC=C2)C")
    ]
    for name, smi in test_smiles:
        result, reason = is_2_5_diketopiperazines(smi)
        print(f"{name}: {result}\n  Reason: {reason}")
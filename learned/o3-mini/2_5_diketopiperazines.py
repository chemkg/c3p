"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:2,5-diketopiperazines
Definition: Any piperazinone that has a piperazine-2,5-dione skeleton.
The following implementation uses an improved SMARTS approach.
The core ring should consist of two nitrogen atoms and four carbon atoms,
with two of the carbon atoms (at positions corresponding to 2 and 5) bearing a carbonyl.
The SMARTS pattern is defined as:
    [#7;R]1[#6;R](=O)[#6;R][#7;R][#6;R](=O)[#6;R]1
which allows additional substituents or stereochemical annotations without preventing a match.
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) contains a 2,5-diketopiperazine core.
    That is, it looks for a six-membered ring which contains two nitrogen atoms and two carbonyl groups 
    (carbon atoms that are double-bonded to oxygen), corresponding to the piperazine-2,5-dione skeleton.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a piperazine-2,5-dione skeleton.
        str: A reason for the classification decision.
    """

    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define a SMARTS pattern for the 2,5-diketopiperazine core.
    # The pattern selects a 6-membered ring with:
    #   - Two nitrogen atoms (atomic number 7) in the ring.
    #   - Four carbon atoms (atomic number 6) in the ring,
    #     two of which are carbonyl carbons (i.e. are double-bonded to oxygen).
    diketopiperazine_smarts = "[#7;R]1[#6;R](=O)[#6;R][#7;R][#6;R](=O)[#6;R]1"
    pattern = Chem.MolFromSmarts(diketopiperazine_smarts)
    if pattern is None:
        return False, "Error constructing the SMARTS pattern."

    # Search for the pattern in the molecule, ignore chirality for a broader match.
    matches = mol.GetSubstructMatches(pattern, useChirality=False)
    if matches:
        # Check that at least one match consists of six atoms (i.e. the complete core ring)
        for match in matches:
            if len(match) == 6:
                return True, "Molecule contains a piperazine-2,5-dione skeleton."
        return False, "Substructure found, but does not form a complete 6-membered diketopiperazine ring."
    else:
        return False, "Molecule does not contain a piperazine-2,5-dione skeleton."

# Example usage (you can add more tests as needed):
if __name__ == "__main__":
    test_molecules = [
        ("Brocazine F", "S1S[C@]23N([C@@H]4[C@@H](O)C=C[C@@H]([C@H]4C2)O)C([C@]15N([C@@H]6[C@@H](O)C=CC([C@H]6C5)=O)C3=O)=O"),
        ("mycocyclosin", "Oc1ccc2C[C@@H]3NC(=O)[C@H](Cc4ccc(O)c(c4)-c1c2)NC3=O"),
        ("tardioxopiperazine A", "O=C1N[C@H](C(=O)N[C@H]1CC=2C=3C(NC2C(C=C)(C)C)=CC=C(C3)CC=C(C)C)C"),
        ("piperazine-2,5-dione", "O=C1CNC(=O)CN1")
    ]
    
    for name, smi in test_molecules:
        result, reason = is_2_5_diketopiperazines(smi)
        print(f"{name}: {result}\n  Reason: {reason}")
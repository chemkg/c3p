"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: Aliphatic Nitrile
Definition: Any nitrile derived from an aliphatic compound.
This program checks for the presence of a nitrile (C#N) group and then
verifies that the nitrile carbon is not aromatic and is attached to a non‐aromatic (aliphatic) group.
"""

from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile must contain at least one nitrile (C#N) group,
    in which the nitrile carbon is not aromatic (and not attached to an aromatic system).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise.
        str: Reason for the classification.
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a nitrile group: a carbon triple bonded to a nitrogen.
    # We first use a simple pattern "[#6]#[#7]" to find nitrile groups.
    nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group found in the molecule"

    # Check each nitrile match
    for match in nitrile_matches:
        # The match returns a tuple of atom indices, where the first is the carbon and the second is nitrogen.
        nitrile_carbon = mol.GetAtomWithIdx(match[0])
        nitrile_nitrogen = mol.GetAtomWithIdx(match[1])

        # Check if the nitrile carbon is aromatic. If yes, then it is not aliphatic.
        if nitrile_carbon.GetIsAromatic():
            return False, "Nitrile carbon is aromatic"

        # Get neighbors of the nitrile carbon (other than the nitrile nitrogen)
        # In a typical nitrile, the nitrile carbon has one other substituent (the R group),
        # so we check that substituent to ensure it is not aromatic.
        for neighbor in nitrile_carbon.GetNeighbors():
            if neighbor.GetIdx() == nitrile_nitrogen.GetIdx():
                continue
            if neighbor.GetIsAromatic():
                return False, "Nitrile group is attached to an aromatic system"
        # You could add more checks here if needed (e.g., verifying the type of attachment) 

    # If at least one nitrile group is present and passed all the checks, we conclude
    # that the compound is an aliphatic nitrile.
    return True, "Contains nitrile group with aliphatic environment"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CC#N",  # acetonitrile (should be aliphatic)
        "c1ccccc1C#N",  # benzonitrile (aromatic nitrile, should fail)
        "CC(C)CC#N"  # isovaleronitrile (aliphatic)
    ]
    for smi in test_smiles:
        result, reason = is_aliphatic_nitrile(smi)
        print(f"SMILES: {smi} -> {result}, {reason}")
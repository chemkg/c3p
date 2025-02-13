"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
# aliphatic_aldoxime.py
"""
Classifies: Any aldoxime derived from an aliphatic aldehyde.
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is any aldoxime (C=N–OH) that is derived from an aliphatic aldehyde.
    This means that the carbon of the aldoxime (originally the aldehyde carbon) must be a CH unit,
    indicating its derivation from an aldehyde (and not a ketone), and should not be part of an aromatic system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic aldoxime, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to ensure that CH units are correctly recognized
    mol = Chem.AddHs(mol)

    # Define a SMARTS pattern for the aldoxime group.
    # This pattern captures a non-aromatic carbon with exactly one hydrogen (CH1) double bonded to a nitrogen
    # that is connected to an oxygen with an explicit H.
    aldoxime_smarts = "[CH1;!a]=[N][O;H]"
    aldoxime_pattern = Chem.MolFromSmarts(aldoxime_smarts)
    if aldoxime_pattern is None:
        return False, "Error in SMARTS pattern"

    # Search for the aldoxime substructure in the molecule.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime group (C=N–OH with a CH unit) found"

    # For each match, check that the origin aldehyde carbon is aliphatic.
    # In an aldehyde, the carbon attached to the C=N–OH group should either have no other heavy atom neighbor
    # or just one that is non-aromatic (or it could be formaldehyde where the only substituents are H).
    for match in matches:
        # match: (carbon_index, nitrogen_index, oxygen_index)
        carbon_idx = match[0]
        nitrogen_idx = match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # Ensure the carbon is not in an aromatic system
        if carbon_atom.GetIsAromatic():
            continue  # this aldoxime appears to come from an aromatic system

        # Now check the neighbors (excluding the nitrogen which is part of the oxime group)
        is_aliphatic = True
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIdx() == nitrogen_idx:
                continue
            # If the neighbor is a heavy atom (not hydrogen) and aromatic, then the originating aldehyde is aromatic.
            if neighbor.GetAtomicNum() != 1 and neighbor.GetIsAromatic():
                is_aliphatic = False
                break

        if is_aliphatic:
            return True, "Aldoxime group derived from an aliphatic aldehyde found"

    return False, "Aldoxime group found, but it appears to be derived from an aromatic aldehyde"
  
# Optional: For testing purposes, you can call the function with one of the example SMILES strings.
if __name__ == "__main__":
    test_smiles = "[H]\\C(C)=N\\O"  # (Z)-acetaldehyde oxime
    result, reason = is_aliphatic_aldoxime(test_smiles)
    print(result, "->", reason)
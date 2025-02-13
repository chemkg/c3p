"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
A compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.
"""

from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule contains an amine functional group.
    An amine is defined as a compound derived from ammonia (NH3) by replacing one, two or three hydrogen atoms by hydrocarbyl groups.
    This function identifies a nitrogen atom whose total substituent count (neighbors + implicit/exlicit hydrogens) equals 3.
    It intentionally excludes nitrogen atoms that are part of an amide (N attached to a carbonyl group)
    or are quaternary (having four substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if at least one amine group is found, False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms to search for relevant nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:   # skip non-nitrogen atoms
            continue

        # Exclude nitrogen with an explicit formal charge that might indicate a quaternary ammonium
        # (or other nitrogen not derived directly from ammonia by substitution)
        if atom.GetFormalCharge() > 0:
            continue

        # Calculate the total number of substituents derived from bonds and implicit hydrogens.
        # For example, in a tertiary amine, there are three heavy atom neighbors and 0 hydrogens.
        total_substituents = atom.GetTotalNumHs() + len(atom.GetNeighbors())

        # The amine definition expects the nitrogen to have three substituents (as in NH3, RNH2, R2NH, or R3N).
        if total_substituents != 3:
            continue

        # Exclude if the nitrogen is directly bonded to a carbon that is part of a carbonyl group (i.e. amide formation)
        amide_flag = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # carbon neighbor
                # Look at bonds from the neighbor to oxygen; if any of these is a double bond then exclude.
                for bond in neighbor.GetBonds():
                    # Check if the bond is a double bond and if the other atom is oxygen.
                    if bond.GetBondTypeAsDouble() and bond.GetOtherAtom(neighbor).GetAtomicNum() == 8:
                        amide_flag = True
                        break
            if amide_flag:
                break

        if amide_flag:
            continue

        # If we have reached here then we have a nitrogen with exactly 3 substituents and no carbonyl involvement.
        return True, "Found at least one amine functional group derived from ammonia"

    # If no such nitrogen is found then the molecule does not conform to the amine definition.
    return False, "No amine functional group (derived from ammonia) identified"
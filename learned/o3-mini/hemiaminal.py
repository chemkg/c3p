"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemiaminal
Definition: “Any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom.
Hemiaminals are intermediates in the formation of imines by addition of an amine to an aldehyde or ketone; 
those derived from primary amines are particularly unstable.”
 
This version explicitly adds hydrogens and then iterates through each saturated (sp3) carbon that “should”
have four substituents. For each such carbon the algorithm checks whether:
  – It has a hydroxyl group (O with an explicit hydrogen, and the bond is a single bond)
  – It has an amino group (N that is not part of an amide—a simple check is to ensure the N is not 
    double bonded to an oxygen)
If found, the function returns True along with the matching atom indices.
"""

from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has a tetrahedral (mostly sp3) carbon that bears both a hydroxyl group (-OH)
    and an amino group (-NH2, -NHR, or -NR2) attached directly to that same carbon.

    The approach is:
      1. Parse the molecule and add explicit hydrogens.
      2. Iterate over sp3 carbon atoms that have a total of 4 substituents.
      3. For each such carbon, check if one neighbor is an oxygen that is in an -OH group
         (a single bond and at least one explicit hydrogen), and at least one neighbor is a
         nitrogen whose bonds do not indicate an amide character (i.e. no double bond to O).
      4. If any carbon satisfies these conditions, return True along with the match details.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a hemiaminal is detected, False otherwise.
        str: Reason/message with matched atom indices or explanation.
    """
    # Parse SMILES. If invalid, return appropriate message.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string; could not parse molecule"

    # Add explicit hydrogens so we can test for hydroxyls properly.
    mol = Chem.AddHs(mol)

    hemiaminal_matches = []  # we will collect tuples (carbon_index, oxygen_index, nitrogen_index)

    # Iterate over atoms to find candidate carbon centers.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # look only at carbons
            continue
        # Check if the carbon is sp3 and has four substituents in total (including implicit H,
        # but since we added Hs all should now be explicit).
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        if atom.GetDegree() != 4:
            continue

        # For each candidate carbon, check its neighbors.
        neighbors = atom.GetNeighbors()
        o_matches = []  # oxygen neighbors representing possible -OH groups
        n_matches = []  # nitrogen neighbors representing possible amino groups
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # We require a single (sigma) bond
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue

            # Look for an oxygen that is an -OH (we require at least one hydrogen on it)
            if nbr.GetAtomicNum() == 8:
                # Check that the oxygen has at least one hydrogen attached.
                # (Note: after AddHs, hydrogens are explicit)
                has_hydrogen = any(neigh.GetAtomicNum() == 1 for neigh in nbr.GetNeighbors())
                if has_hydrogen:
                    o_matches.append(nbr.GetIdx())
            # Look for a nitrogen and try to rule out amide-like nitrogens.
            elif nbr.GetAtomicNum() == 7:
                # For the nitrogen neighbor, ensure it does not have any double bond to oxygen
                # (A double bond to oxygen might indicate an amide or similar functionality).
                is_amide = False
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetAtomicNum() == 8:
                        b = mol.GetBondBetweenAtoms(nbr.GetIdx(), nnbr.GetIdx())
                        if b is not None and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            is_amide = True
                            break
                if not is_amide:
                    n_matches.append(nbr.GetIdx())

        # If we found at least one OH group and one amino group, record them.
        if o_matches and n_matches:
            # Record all combinations from this carbon for clarity.
            for o_idx in o_matches:
                for n_idx in n_matches:
                    hemiaminal_matches.append((atom.GetIdx(), o_idx, n_idx))

    if hemiaminal_matches:
        return True, f"Found hemiaminal substructure in atoms with indices: {hemiaminal_matches}"
    else:
        return False, "No carbon with both hydroxyl and amino substituents found; not a hemiaminal"

# Example usage
if __name__ == "__main__":
    # You may test with one of the provided SMILES strings.
    test_smiles = "OC(N)CC"  # 2-Aminopropanol should be a simple hemiaminal
    result, reason = is_hemiaminal(test_smiles)
    print(result, reason)
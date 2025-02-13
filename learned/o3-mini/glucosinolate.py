"""
Classifies: CHEBI:24279 glucosinolate
"""
#!/usr/bin/env python
"""
Classifies: Glucosinolate

Definition:
  Glucosinolates are water‐soluble anionic substituted thioglucosides.
  They have a central C atom bonded via an S atom to a glycone (sugar) group,
  and via an N atom to a sulfonated oxime group, in addition to carrying a side‐group.
  The side-chain and sulfate group are in an anti (E) configuration across the C=N double bond.
  
This program tries to detect that pattern by:
  - Looking for a pyranose (glucose-like) sugar motif.
  - Searching for a central motif defined by a sulfur connected to a C which is double-bonded to an N,
    where that N is attached to an oxygen that connects to a sulfonated sulfate group.
  - Verifying that the central carbon also carries an extra substituent (side chain).
  - Checking that the C=N double bond has the explicit E (anti) stereochemistry.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondStereo

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.

    The algorithm checks that:
      - There is a glycone (sugar) moiety.
      - There exists a central motif: an S atom attached to a central C that is double-bonded to an N,
        with that N connected to an O which in turn is bonded to a sulfonate group (S(=O)(=O)[O-]).
      - The central C bears an extra substituent (side chain).
      - The C=N double bond stereochemistry is explicitly E (anti).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule meets glucosinolate criteria, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the glycone (sugar) moiety.
    # Use a general pyranose pattern without chirality markers.
    sugar_smarts = "C1OC(C(O)C(O)C1O)CO"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "Glycone (sugar) moiety not found"

    # Define a SMARTS for the central glucosinolate motif.
    # This pattern describes:
    #  - A tetrahedral S atom (from the thiol linkage) connected to a central C,
    #  - The central C is double bonded to an N,
    #  - The N is bound to an O which is attached to a sulfonate group.
    # Note: We do not encode the side-chain in the SMARTS; we instead check that the central C has an extra substituent.
    central_smarts = "[S;X2]-[C](=[N]-[O]-[S](=O)(=O)[O-])"
    central_pattern = Chem.MolFromSmarts(central_smarts)
    central_matches = mol.GetSubstructMatches(central_pattern)
    if not central_matches:
        return False, "Central glucosinolate motif not found"

    # Loop over possible matches for the central motif.
    valid = False
    for match in central_matches:
        # Expected match order: [S, C, N, O, S(sulfate)]
        s_index = match[0]
        c_index = match[1]
        n_index = match[2]
        # (Indices 3 and 4 are part of the oxime-sulfonate group.)

        # Check that the S atom (thioglucoside) is attached to a glycone (sugar) moiety.
        s_atom = mol.GetAtomWithIdx(s_index)
        sugar_connected = False
        for nbr in s_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            for sugar in sugar_matches:
                if nbr_idx in sugar:
                    sugar_connected = True
                    break
            if sugar_connected:
                break
        if not sugar_connected:
            continue  # This motif is not attached to a sugar.

        # Verify that the central carbon has at least three neighbors;
        # it should be bonded to the S atom, the N atom (double bond), and a side-chain.
        c_atom = mol.GetAtomWithIdx(c_index)
        if c_atom.GetDegree() < 3:
            continue

        # Check the stereochemistry of the C=N double bond.
        bond = mol.GetBondBetweenAtoms(c_index, n_index)
        if bond is None:
            continue
        stereo = bond.GetStereo()
        if stereo != BondStereo.STEREOE:
            continue

        # If all tests pass for this match, classify the molecule as a glucosinolate.
        valid = True
        break

    if valid:
        return True, "Molecule contains glucosinolate motif with glycone connection and anti C=N configuration"
    else:
        return False, "Central motif found but either not attached to a glycone, side-chain missing, or the C=N bond is not E configured"

# For testing purposes (optional)
if __name__ == "__main__":
    test_smiles = [
        "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CCC=C",  # Example 1
        "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS([O-])(=O)=O)/CCCCC",  # pentylglucosinolate
        "S(C1OC(C(O)C(O)C1O)CO)/C(=N\\OS([O-])(=O)=O)/CCCC",  # butylglucosinolate: likely to test stereochemistry & sugar attachment
    ]
    
    for s in test_smiles:
        result, reason = is_glucosinolate(s)
        print("SMILES:", s)
        print("Glucosinolate:", result, "|", reason)
        print("----")
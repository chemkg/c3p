"""
Classifies: CHEBI:24279 glucosinolate
"""
#!/usr/bin/env python
"""
Classifies: Glucosinolate
Definition:
  Glucosinolates are water‐soluble anionic substituted thioglucosides.
  They have a central C atom which is bonded via an S atom to a glycone (sugar) group,
  and via an N atom to a sulfonated oxime group, and additionally carry a side‐group.
  The side-chain and sulfate group are in an anti (E) configuration across the C=N double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondStereo

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    The algorithm searches for:
      - A glycone (sugar) moiety. Here we use a generalized pyranose SMARTS pattern,
        chosen to match many glucose-like rings: "C1OC(C(O)C(O)C1O)CO".
      - A central motif consisting of an S atom bound to a C that is double bonded
        to an N, which in turn is attached to an O that bears a sulfonate group.
      - An explicit anti (E) stereochemical configuration for the C=N double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a glucosinolate, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a more general SMARTS for a pyranose (glucose-like) ring.
    # This pattern is chosen to capture many common ways the sugar moiety is written.
    sugar_smarts = "C1OC(C(O)C(O)C1O)CO"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "Glycone (sugar) moiety not found"
        
    # Define a SMARTS for the glucosinolate central motif:
    # [S]-C(=N-O-[S](=O)(=O)[O-])
    # Note: this does not encode stereochemistry; that is checked explicitly later.
    central_smarts = "[S]-C(=N-O-[S](=O)(=O)[O-])"
    central_pattern = Chem.MolFromSmarts(central_smarts)
    central_matches = mol.GetSubstructMatches(central_pattern)
    if not central_matches:
        return False, "Central glucosinolate motif not found"
        
    # Loop over each match of the central motif.
    # The expected match tuple is of the order: (i_S, i_C, i_N, i_O, i_Sulfate).
    found_valid = False
    for match in central_matches:
        s_index = match[0]
        c_index = match[1]
        n_index = match[2]
        # Verify that the S atom connecting the central C is attached to a glycone ring.
        s_atom = mol.GetAtomWithIdx(s_index)
        sugar_connected = False
        for nbr in s_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Check if the neighbor belongs to any sugar match.
            for sugar in sugar_matches:
                if nbr_idx in sugar:
                    sugar_connected = True
                    break
            if sugar_connected:
                break
        if not sugar_connected:
            # Skip this match if the glycone is not attached.
            continue

        # Check the stereochemistry of the C=N double bond.
        bond = mol.GetBondBetweenAtoms(c_index, n_index)
        if bond is None:
            continue

        stereo = bond.GetStereo()
        # We require the double bond to have an explicit anti (E) configuration.
        if stereo != BondStereo.STEREOE:
            # If stereo is undefined or not E, then skip this match.
            continue
        
        # If both tests pass for this match, then the molecule is classified as a glucosinolate.
        found_valid = True
        break
        
    if found_valid:
        return True, "Molecule contains glucosinolate motif with glycone connection and anti C=N configuration"
    else:
        return False, "Central motif found but either not attached to a glycone or the C=N bond is not E configured"

# For testing (optional)
if __name__ == "__main__":
    # List some example SMILES of glucosinolates.
    test_smiles = [
        "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CCC=C",  # Example with sugar missing in previous check
        "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS([O-])(=O)=O)/CCCCC",  # pentylglucosinolate
        "S(C1OC(C(O)C(O)C1O)CO)/C(=N\\OS([O-])(=O)=O)/CCCC",  # butylglucosinolate
    ]
    
    for s in test_smiles:
        result, reason = is_glucosinolate(s)
        print("SMILES:", s)
        print("Glucosinolate:", result, "|", reason)
        print("----")
"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: Glucosinolate

Definition:
  Glucosinolates are water‐soluble anionic substituted thioglucosides.
  They have a central C atom bonded via an S atom to a glycone (sugar) group,
  and via an N atom to a sulfonated oxime group, in addition to carrying a side‐group.
  The side-chain and sulfate group are in an anti (E) configuration across the C=N double bond.
  
The improved algorithm:
  - Looks for a generic glycone/sugar moiety (pyranose ring) using a standard glucose pattern.
  - Searches for a central motif defined as an S attached to a C that is double-bonded to an N;
    the N is attached to an O that in turn is connected to a sulfonate group.
    We use a flexible SMARTS for that motif that accepts the sulfate terminus as either [O-] or O.
  - Checks that the central C has a third substituent apart from S and the double-bonded N.
  - Verifies that the C=N double bond is explicitly assigned the E (anti) stereochemistry.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondStereo

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    This function verifies:
      - Presence of a glycone (sugar) moiety, using a generic glucose SMARTS.
      - Existence of a central motif: an S atom bonded to a C that is double-bonded to an N,
        where the N is attached to an O that connects to a sulfonate group. The SMARTS is
        written to allow the terminal oxygen to be either charged or uncharged.
      - That the central C has an extra substituent (side-chain).
      - That the C=N double bond is explicitly E (anti) configured.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets glucosinolate criteria, False otherwise.
        str: Reason for classification decision.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the glycone (pyranose) moiety.
    # This is a generic glucose-like pattern (without requiring exact chirality).
    sugar_smarts = "C1OC(C(O)C(O)C1O)CO"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "Glycone (sugar) moiety not found"

    # Define a flexible SMARTS for the central glucosinolate motif.
    # The pattern encodes:
    #  - A sulfur atom (S) 
    #  - Connected to a carbon (C) that is double-bonded to a nitrogen (N)
    #  - The nitrogen is single bonded to an oxygen (O) that is attached to a sulfonate group,
    #    with the terminal oxygen allowed to be either O or [O-].
    # Using atomic numbers to be a bit more generic.
    central_smarts = "[#16]-[#6](=[#7]-[#8]-[#16](=O)(=O)[O,$([O-])])"
    central_pattern = Chem.MolFromSmarts(central_smarts)
    central_matches = mol.GetSubstructMatches(central_pattern)
    if not central_matches:
        return False, "Central glucosinolate motif not found"

    # Loop through all central motif matches to check for sugar connection, extra substituent,
    # and explicit E (anti) stereochemistry of the C=N double bond.
    for match in central_matches:
        # Expected match order according to the SMARTS: 
        # index0 -> S, index1 -> C (central), index2 -> N, index3 -> O, index4 -> S of sulfonate.
        s_idx = match[0]
        c_idx = match[1]
        n_idx = match[2]
        
        # Check that the S atom is connected to the glycone (sugar) moiety.
        s_atom = mol.GetAtomWithIdx(s_idx)
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
            # This central motif is not linked to a sugar; try next match.
            continue

        # Check that the central carbon has at least three neighbors 
        # (S, double-bonded N, and an extra side chain)
        c_atom = mol.GetAtomWithIdx(c_idx)
        if c_atom.GetDegree() < 3:
            continue

        # Check the stereochemistry of the C=N double bond.
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond is None:
            continue
        stereo = bond.GetStereo()
        if stereo != BondStereo.STEREOE:
            continue

        # If all conditions are met, return positive classification.
        return True, "Molecule contains glucosinolate motif with glycone connection and anti C=N configuration"
    
    # If we exit the loop without a valid match:
    return False, "Central motif found but either not attached to a glycone, side-chain missing, or the C=N bond is not E configured"

# For optional testing
if __name__ == "__main__":
    test_smiles = [
        "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CCC=C",  # {(e)-...} glucosinolate
        "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS([O-])(=O)=O)/CCCCC",  # pentylglucosinolate
        "S(C1OC(C(O)C(O)C1O)CO)/C(=N\\OS([O-])(=O)=O)/CCCC",  # butylglucosinolate (stereochemistry may be absent/inexact)
    ]
    
    for s in test_smiles:
        result, reason = is_glucosinolate(s)
        print("SMILES:", s)
        print("Glucosinolate:", result, "|", reason)
        print("----")
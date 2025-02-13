"""
Classifies: CHEBI:24279 glucosinolate
"""
#!/usr/bin/env python
"""
Classifies: Glucosinolate
Definition: Water‐soluble anionic substituted thioglucosides.
Glucosinolates have a central C atom which is bonded via an S atom to a glycone group 
and via an N atom to a sulfonated oxime group – and which also carries a side‐group.
The side‐chain and the sulfate group have an anti stereochemical configuration 
(across the C=N double bond, i.e. “E” configuration).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondStereo

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    The algorithm checks for the presence of a pyranose sugar (glycone)
    and the central motif "S-C(=N-O-S(=O)(=O)[O-])". In this motif the S atom is 
    connected to the sugar moiety and the central C that is double bonded to an N 
    that in turn is attached to an OSO3– group. Finally, it verifies that the C=N double bond 
    exhibits an anti (E) stereochemical configuration.
    
    Note: The approach is based on substructure SMARTS matching and relies on 
    explicit stereochemical information in the SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a glucosinolate, False otherwise
        str: Reason for the classification result
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for a pyranose (six-membered sugar) ring typical of glucose.
    # This pattern ignores exact stereochemistry to allow for minor variations.
    sugar_smarts = "C1OC(O)C(O)C(O)C(O)C1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "Glycone (sugar) moiety not found"
    
    # Define a SMARTS for the central oxime motif.
    # The pattern is written as:
    # [S]-C(=N-O-[S](=O)(=O)[O-])
    # That is, an S atom attached to a C that is double-bonded to an N, 
    # which in turn is single-bonded to an O that is connected to a sulfonate group.
    central_smarts = "[S]-C(=N-O-[S](=O)(=O)[O-])"
    central_pattern = Chem.MolFromSmarts(central_smarts)
    central_matches = mol.GetSubstructMatches(central_pattern)
    if not central_matches:
        return False, "Central glucosinolate motif not found"

    # For each match in the central motif, verify:
    # (a) the S (first atom in the match) is connected to a sugar ring,
    # (b) the C=N bond has anti (E) stereochemistry.
    for match in central_matches:
        # match is a tuple of atom indices corresponding to:
        # (i_S, i_C, i_N, i_O, i_S_sulfate)
        s_index = match[0]
        c_index = match[1]
        n_index = match[2]
        # Check that the S atom (the one connecting to the central C) is bound to a sugar.
        s_atom = mol.GetAtomWithIdx(s_index)
        sugar_connected = False
        for nbr in s_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # See if the neighbor atom index is part of any sugar match.
            for sugar in sugar_matches:
                if nbr_idx in sugar:
                    sugar_connected = True
                    break
            if sugar_connected:
                break
        if not sugar_connected:
            continue  # Try next central match
        
        # Check the stereochemistry of the C=N double bond.
        bond = mol.GetBondBetweenAtoms(c_index, n_index)
        # If stereo is not defined, we treat that as a failure.
        if bond is None:
            continue
        
        stereo = bond.GetStereo()
        if stereo != BondStereo.STEREOE:
            # Sometimes the stereo may be undefined (STEREONONE) so we also reject that.
            continue
        
        # If we reached here, we found a central motif attached to a sugar 
        # and the C=N bond is E configured.
        return True, "Molecule contains glucosinolate motif with correct glycone link and anti C=N configuration"
    
    # If none of the central motifs pass both tests:
    if central_matches:
        return False, "Central motif found but either not connected to a sugar moiety or C=N bond is not anti configured"
    
    return False, "Unknown error in glucosinolate classification"

# For testing (optional)
if __name__ == "__main__":
    # Example SMILES for pentylglucosinolate:
    smiles_examples = [
        "S(C1OC(C(O)C(O)C1O)CO)/C(=N\\OS([O-])(=O)=O)/CCCCC",
        "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)SC(=NOS([O-])(=O)=O)C[C@H](C=C)O"
    ]
    for s in smiles_examples:
        result, reason = is_glucosinolate(s)
        print("SMILES:", s)
        print("Glucosinolate:", result, "|", reason)
        print("----")
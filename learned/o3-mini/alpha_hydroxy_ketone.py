"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha‐carbon relative to the C=O group.
Improved approach:
  • Iterate over atoms to find a carbonyl carbon (i.e. a carbon with a double-bonded oxygen)
  • Check that the carbonyl carbon has exactly two other substituents, both being carbon (i.e. a ketone not an aldehyde or carboxylic acid)
  • For each of these neighboring (alpha) carbons, check that it is sp3-hybridized and has at least one neighboring oxygen that bears at least one hydrogen (i.e. an –OH group)
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (C(=O) with two carbon substituents) that has an -OH
    group attached on at least one of the carbons adjacent (alpha-carbon) to the C=O.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an alpha-hydroxy ketone, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms of the molecule to find potential carbonyl carbons.
    for atom in mol.GetAtoms():
        # We are looking for carbon atoms only.
        if atom.GetAtomicNum() != 6:
            continue

        # Look for double bonds to oxygen.
        dbl_oxygen_neighbors = []
        for bond in atom.GetBonds():
            # Check for a double bond to an oxygen.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    dbl_oxygen_neighbors.append(nbr)

        # To be a carbonyl, there must be at least one (typically one) double-bonded oxygen.
        if len(dbl_oxygen_neighbors) < 1:
            continue

        # Now, get the other substituents.
        # For a ketone, apart from the double-bonded oxygen, the carbonyl carbon should have exactly two more groups,
        # and both should be carbon atoms.
        other_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr not in dbl_oxygen_neighbors]
        if len(other_neighbors) != 2:
            continue
        if not all(nbr.GetAtomicNum() == 6 for nbr in other_neighbors):
            continue  # Skip if any substituent is not a carbon.

        # For each alpha carbon (i.e. each of the two substituents),
        # check if it is sp3-hybridized and if it bears an -OH group.
        for alpha in other_neighbors:
            if alpha.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            
            # Look for an oxygen neighbor attached to this alpha carbon that carries at least one hydrogen.
            for neighbor in alpha.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                    return True, ("Molecule contains a ketone (C=O) with an sp3 alpha-carbon bearing a hydroxy group "
                                  "(alpha-hydroxy ketone detected).")
                    
    return False, "No appropriate alpha-hydroxy ketone motif found in molecule."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "OCC(O)C(=O)CO",  # erythrulose, expected True
        "C[C@H]1[C@H]2[C@H](Cc3c[nH]c4ccccc34)NC(=O)[C@@]22[C@@H](\\C=C\\C[C@H](C)\\C=C(C)\\C(=O)[C@@H](O)CCC2=O)[C@@H]2O[C@]12C",  # chaetoglobosin F, complex example
        "(R)-3-hydroxy-3-methyl-2-oxopentanoic acid",  # should be True if the ketone is correctly identified
    ]
    for smi in test_smiles:
        res, reason = is_alpha_hydroxy_ketone(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")
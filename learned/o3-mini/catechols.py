"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: Catechols (compounds containing an o-diphenol component)
A catechol is any compound that contains an aromatic ring with two adjacent hydroxyl groups.
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol requires at least one aromatic ring in the structure where two adjacent
    carbon atoms both have a hydroxyl (-OH) substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an o-diphenol (catechol) component, False otherwise
        str: Reason for classification
    """
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper function: Check if an atom has a hydroxyl substituent.
    # We consider an -OH group when an oxygen atom (atomic number 8) is attached
    # to the atom and has at least one hydrogen (implicit or explicit).
    def has_hydroxyl(atom):
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() > 0:
                # Also check that the oxygen is not doubly bonded to the atom
                # (this avoids misidentifying carbonyl oxygens as -OH)
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() < 2.0:
                    return True
        return False

    # Loop over aromatic carbon atoms: only consider atoms that are aromatic carbons
    # which may serve as the backbone of an aromatic ring.
    for atom in mol.GetAtoms():
        # Check if this atom is aromatic and a carbon and has a hydroxyl group attached
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and has_hydroxyl(atom):
            # Check neighbors: look for an adjacent aromatic carbon also bearing a hydroxyl group.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic() and has_hydroxyl(nbr):
                    # Found two adjacent aromatic carbons with -OH substituents.
                    return True, "Found neighboring aromatic carbons with hydroxyl groups (o-diphenol component)"
    
    return False, "No adjacent hydroxyl groups on any aromatic ring found"

# Example tests (uncomment to run)
# test_smiles = [
#     "OC1=C(O)C=CC=C1CCCC/C=C\\C/C=C\\CCCCCCCC2=C(O)C(O)=CC=C2",  # example structure with catechol component
#     "COc1cc(CCc2ccc(O)c(O)c2)cc(O)c1O",  # dendrocandin E, contains catechol fragment
#     "Nc1ccc(O)c(O)c1",  # 4-Aminocatechol, contains catechol fragment
#     "C1=CC=CC=C1",  # benzene, no hydroxyl groups
# ]
# for sm in test_smiles:
#     flag, reason = is_catechols(sm)
#     print(f"SMILES: {sm}\nClassification: {flag}\nReason: {reason}\n")
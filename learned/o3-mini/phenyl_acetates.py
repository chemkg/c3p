"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
with the hydroxy group of any phenol.
"""

from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester (–OC(=O)CH3) where the ester oxygen is attached to an aromatic (phenol)
    ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the phenyl acetate criteria, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an acetate ester group: O-C(=O)-C (explicitly CH3 in the acetyl group)
    acetate_smarts = "O-C(=O)C"
    acetate_mol = Chem.MolFromSmarts(acetate_smarts)
    if acetate_mol is None:
        return False, "Error in acetate SMARTS pattern"

    # Look for the acetate ester substructure within the molecule
    matches = mol.GetSubstructMatches(acetate_mol)
    if not matches:
        return False, "No acetate ester group found"

    # Check that for at least one match, the oxygen atom is attached to an aromatic ring (phenol part)
    for match in matches:
        # match is a tuple of atom indices corresponding to the SMARTS "O-C(=O)C"
        # In our pattern:
        #   index 0: oxygen of the ester group
        #   index 1: carbonyl carbon
        #   index 2: methyl group (CH3)
        ox_idx = match[0]
        ox_atom = mol.GetAtomWithIdx(ox_idx)
        aromatic_attached = False
        # Loop over neighbors of the oxygen atom
        for nbr in ox_atom.GetNeighbors():
            # We skip the neighbor that belongs to the ester group (carbonyl carbon)
            if nbr.GetIdx() == match[1]:
                continue
            # Check if the neighbor is aromatic and is a carbon (atomic number 6)
            if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                aromatic_attached = True
                break
        if aromatic_attached:
            return True, "Found acetate ester group attached to an aromatic (phenol) ring"
    return False, "Acetate ester group found but not attached to an aromatic ring"
"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Primary Alcohol
Definition: Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon 
          which is itself bonded to an aromatic ring.
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    A qualifying structure contains a primary alcohol group (CH2OH) where the CH2 is 
    directly bound to an aromatic carbon.

    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is an aromatic primary alcohol, False otherwise.
        str: Explanation/reason for the decision.
    """
    # Parse the SMILES string to create a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (ensures accurate hydrogen counts)
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms looking for hydroxyl groups (-OH)
    # We look for an oxygen that is in an –OH group.
    for atom in mol.GetAtoms():
        # Check if the atom is oxygen
        if atom.GetAtomicNum() != 8:
            continue
        
        # Get all neighbors of the oxygen
        neighbors = atom.GetNeighbors()
        
        # A hydroxyl oxygen in –OH (not in carboxylic acids etc.) should have one heavy atom (carbon)
        # and at least one hydrogen.
        heavy_neighbors = [n for n in neighbors if n.GetAtomicNum() > 1]
        hydrogen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 1]
        
        # We require exactly one heavy neighbor (the carbon) and at least one hydrogen
        if len(heavy_neighbors) != 1 or len(hydrogen_neighbors) < 1:
            continue  # not a standard -OH group
        
        # The heavy neighbor must be a carbon.
        carbon = heavy_neighbors[0]
        if carbon.GetAtomicNum() != 6:
            continue
        
        # For a primary alcohol, the carbon should be CH2.
        # Count heavy atom neighbors of the carbon excluding the oxygen we came from.
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() > 1 and n.GetIdx() != atom.GetIdx()]
        
        # For a primary alcohol (CH2) the carbon has exactly one heavy neighbor aside from the -OH.
        if len(carbon_neighbors) != 1:
            continue  # Not a primary alcohol, might be secondary or tertiary
        
        # Check that the one heavy neighbor (aside from -OH) is part of an aromatic ring.
        # RDKit flags an atom as aromatic via GetIsAromatic().
        other_neighbor = carbon_neighbors[0]
        if not other_neighbor.GetIsAromatic():
            continue  # The carbon is not bonded to an aromatic ring
        
        # If we reached here, we found a valid aromatic primary alcohol group.
        return True, "Found an aromatic primary alcohol group (CH2OH attached to an aromatic ring)"

    # If we looped through all oxygen atoms and found none qualifying:
    return False, "No aromatic primary alcohol group found in the molecule"


# Some example testing (can be removed or commented out when using as a module):
if __name__ == "__main__":
    test_smiles = [
        "CCCC1=C(C(O)=O)C(CO)=CO1",  # 4-(hydroxymethyl)-2-propylfuran-3-carboxylic acid (contains the group? might have extra functionalities though)
        "CC1=CC=C(CO)C(C)=C1",      # 2,4-dimethylbenzyl alcohol - should return True
        "C1=CC(=CN=C1)CO",          # 3-pyridinemethanol - should return True
        "CCCC",                     # simple alkane, should return False
    ]
    for smi in test_smiles:
        result, reason = is_aromatic_primary_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
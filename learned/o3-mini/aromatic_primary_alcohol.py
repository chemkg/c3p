"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Primary Alcohol
Definition: Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon 
            which is itself bonded to an aromatic ring.
In practice, we are looking for an aliphatic CH2-group that forms a –CH2OH moiety where:
  - The carbon is non‐aromatic.
  - It has exactly 2 hydrogens (making it a primary carbon).
  - It is bonded to exactly 2 heavy atoms: one oxygen (which should be in an –OH group) 
    and one aromatic heavy atom.
This implementation iterates over all carbons in the molecule (after adding explicit hydrogens)
and checks if any group meets these precise conditions.
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains an aromatic primary alcohol group.
    A qualifying structure is one in which a CH2 unit (with exactly two hydrogens) is bonded
    to exactly two heavy atoms: one being an -OH oxygen (i.e. an oxygen attached to at least one 
    hydrogen) and the other an aromatic atom. This corresponds to an Ar-CH2OH moiety.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule contains an aromatic primary alcohol group, False otherwise.
        str: Explanation/reason for the decision.
    """
    # Convert the SMILES to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are accurate
    mol = Chem.AddHs(mol)
    
    # Iterate over every carbon atom in the molecule
    for atom in mol.GetAtoms():
        # We are interested only in carbon atoms (atomic number 6) which are not aromatic;
        # in benzyl alcohol, the CH2 group is not part of the aromatic ring.
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            continue
        
        # Check if this carbon is CH2: it must have exactly two hydrogen neighbors.
        hydrogen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
        if len(hydrogen_neighbors) != 2:
            continue
        
        # Get heavy neighbors (non-hydrogen atoms) that are bonded to this carbon.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # For a primary alcohol unit (-CH2OH) there must be exactly two heavy neighbors:
        # one should be the hydroxyl oxygen and one should be an aromatic atom.
        if len(heavy_neighbors) != 2:
            continue
        
        # Identify the hydroxyl oxygen and the aromatic neighbor
        hydroxyl_oxygen = None
        aromatic_neighbor = None
        for nbr in heavy_neighbors:
            if nbr.GetAtomicNum() == 8:
                # Verify that the oxygen is in an -OH group: at least one hydrogen attached.
                o_h_count = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                if o_h_count >= 1:
                    hydroxyl_oxygen = nbr
            else:
                # The other neighbor should be aromatic.
                if nbr.GetIsAromatic():
                    aromatic_neighbor = nbr
        
        # If we found both an -OH oxygen and an aromatic neighbor, that's our moiety.
        if hydroxyl_oxygen is not None and aromatic_neighbor is not None:
            return True, "Found an aromatic primary alcohol group (CH2OH attached directly to an aromatic ring)"
    
    # If no matching group is found in the entire molecule
    return False, "No aromatic primary alcohol group found in the molecule"


# For testing the function (you may remove or modify these tests as needed)
if __name__ == "__main__":
    test_smiles = [
        "CCCC1=C(C(O)=O)C(CO)=CO1",   # 4-(hydroxymethyl)-2-propylfuran-3-carboxylic acid
        "CCCCn1c(CO)ccc1C=O",         # 1-butylpyrraline
        "CC1=CC=C(CO)C(C)=C1",        # 2,4-dimethylbenzyl alcohol
        "C1=CC(=CN=C1)CO",            # 3-pyridinemethanol
        "CCC",                       # simple alkane (should be False)
        "C=1C2=C(C=CC1)SC(=C2)C(=O)N[C@@H](CC(C)C)C(=O)N3CCN(CC3)C(=O)[C@H](CO)NS(C=4C=CC(=CC4Cl)Cl)(=O)=O"  # GSK1016790A (previous false negative)
    ]
    
    for smi in test_smiles:
        result, reason = is_aromatic_primary_alcohol(smi)
        print(f"SMILES: {smi}")
        print(f"Result: {result}")
        print(f"Reason: {reason}\n")
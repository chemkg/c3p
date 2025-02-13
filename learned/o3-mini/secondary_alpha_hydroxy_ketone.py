"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary alpha-hydroxy ketone (acyloin)

A secondary alpha-hydroxy ketone is defined by the motif
  R–C(=O)–CH(OH)–R′
i.e. an acyloin center in which the hydroxy-bearing carbon (CH)
bears exactly one hydrogen and is attached to two carbon groups:
one from a carbonyl (C=O) and one organanyl (R′).
This program parses the molecule from a SMILES string,
adds explicit hydrogens, and then checks every carbon atom
to find such an acyloin motif.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin)
    based on its SMILES string.
    
    The algorithm examines each carbon atom for the following:
      - It is a carbon (atomic number 6).
      - It has exactly one hydrogen attached.
      - It has exactly three heavy (non-hydrogen) neighbors.
      - Exactly one of these neighbors is an oxygen attached as –OH 
        (i.e. the oxygen has at least one hydrogen itself).
      - The other two neighbors are carbons and at least one of them
        is a carbonyl carbon (has a double bond to an oxygen).
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a secondary α–hydroxy ketone (acyloin) motif, False otherwise.
        str: A message explaining the result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the count of hydrogens on atoms is reliable
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms to look for a candidate acyloin center.
    for atom in mol.GetAtoms():
        # Look for carbon atoms (atomic number 6)
        if atom.GetAtomicNum() != 6:
            continue
        
        # Count hydrogen neighbors (explicit H's)
        h_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
        # To be secondary, the carbon must have exactly one hydrogen.
        if h_count != 1:
            continue
        
        # Count heavy (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # For a secondary acyloin center, there should be exactly 3 heavy neighbors:
        # one is the –OH group and the other two are carbons (one being the carbonyl carbon)
        if len(heavy_neighbors) != 3:
            continue
        
        # Among the heavy neighbors, identify the -OH neighbor(s).
        oh_neighbors = []
        for nbr in heavy_neighbors:
            if nbr.GetAtomicNum() == 8:
                # Check if this oxygen is in a hydroxyl group by verifying that it has at least one H attached.
                nbr_h = sum(1 for x in nbr.GetNeighbors() if x.GetAtomicNum() == 1)
                if nbr_h >= 1:
                    oh_neighbors.append(nbr)
        if len(oh_neighbors) != 1:
            continue
        
        # Among the other heavy neighbors, find carbon atoms.
        carbon_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 2:
            continue
        
        # Check that at least one of these carbons is a carbonyl carbon:
        # a carbon that is double bonded to an oxygen.
        carbonyl_found = False
        for c_nbr in carbon_neighbors:
            # Loop over bonds for this neighbor
            for bond in c_nbr.GetBonds():
                # Check if the bond is a double bond and the other atom is oxygen.
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(c_nbr)
                    if other.GetAtomicNum() == 8:
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        
        # If a carbonyl neighbor is found, we have the acyloin motif.
        if carbonyl_found:
            return True, "Contains secondary alpha-hydroxy ketone (acyloin) motif"
    
    return False, "Does not contain the secondary alpha-hydroxy ketone (acyloin) motif"

# Example usage (for testing the function)
if __name__ == "__main__":
    # Example SMILES strings for compounds reported to belong to this class.
    test_smiles = [
        "OC(C(=O)c1ccccc1)c1ccccc1",  # benzoin: expected True
        "OC[C@@H](O)C(=O)CO",          # D-erythrulose: expected True
        "O=C1/C(=C(/O)\\C=C\\C=C\\C)/[C@H]2[C@@](O)(C)C([C@@]1(C(OCC)C2)C)=O"  # Rezishanone C: expected False or borderline
    ]
    for smi in test_smiles:
        result, reason = is_secondary_alpha_hydroxy_ketone(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
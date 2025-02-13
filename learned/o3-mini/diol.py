"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol
A diol is defined here as a compound that contains exactly two alcoholic hydroxy groups.
An alcoholic hydroxy group is an –OH where the oxygen is directly bonded to a non‐aromatic carbon,
and that carbon is not part of a carbonyl (C=O). In addition, the two –OH groups should be found on carbons
that are “close” in the molecular framework (i.e., separated by no more than 10 bonds), to avoid classifying
large multifunctional compounds as diols.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    Our definition: a diol contains exactly 2 alcoholic hydroxy groups (–OH groups not part of
    carboxylic acids or attached to aromatic carbons). In addition, the two alcohol groups must
    reside on carbons that are not too distant in the molecular graph.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is (classified as) a diol, False otherwise.
        str: Reason for the classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, find ALL hydroxyl groups (-OH) using SMARTS.
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxy_pattern is None:
        return False, "Failed to create hydroxy group pattern"
        
    all_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # We'll collect indices of alcoholic hydroxyl groups.
    # An alcoholic hydroxyl group is one where the oxygen is bonded to exactly one heavy (non-H) neighbor,
    # that neighbor is carbon, non‐aromatic, and not double‐bonded to an oxygen (i.e. not part of a carbonyl).
    alcoholic_oh_indices = []
    # Also keep the index of the C to which the -OH is attached
    alcohol_carbon_indices = []
    
    for match in all_matches:
        # match is a tuple containing the index of the O atom
        o_atom = mol.GetAtomWithIdx(match[0])
        # find heavy neighbors (non-hydrogen)
        heavy_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            continue  # unusual – skip if not exactly one heavy neighbor
        c_atom = heavy_neighbors[0]
        # Only consider if neighbor is carbon.
        if c_atom.GetAtomicNum() != 6:
            continue
        # Exclude if the carbon is aromatic
        if c_atom.GetIsAromatic():
            continue
        # Exclude if the carbon is likely part of a carbonyl.
        is_carbonyl = False
        for bond in c_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(c_atom)
                if other.GetAtomicNum() == 8:
                    is_carbonyl = True
                    break
        if is_carbonyl:
            continue
        
        # If passed all filters, count this as an alcoholic -OH.
        alcoholic_oh_indices.append(o_atom.GetIdx())
        alcohol_carbon_indices.append(c_atom.GetIdx())
    
    count_alcoholic = len(alcoholic_oh_indices)
    if count_alcoholic != 2:
        return False, f"Found {count_alcoholic} alcoholic hydroxy group(s); a diol requires exactly 2"
    
    # As an extra check, compute the shortest path (in number of bonds)
    # between the two carbons bearing the alcoholic OH groups.
    # In a typical diol these two groups are not very remote.
    try:
        path = Chem.rdmolops.GetShortestPath(mol, alcohol_carbon_indices[0], alcohol_carbon_indices[1])
    except Exception as e:
        return False, f"Error computing shortest path between alcohol-bearing carbons: {e}"
    
    # Set a threshold for bond length between the two alcohol groups.
    # (For many diols the two -OH groups are on adjacent or nearby carbons.)
    max_bond_distance = 10
    if len(path)-1 > max_bond_distance:
        return False, f"Alcoholic hydroxy groups are too far apart (distance {len(path)-1} bonds; threshold is {max_bond_distance})"
    
    return True, "Contains exactly two alcoholic hydroxy groups (diol) with appropriate connectivity"

# Example usage: (you can remove or update this section if only the function is needed)
if __name__ == "__main__":
    # Test example: (R,R)-butane-2,3-diol
    example_smiles = "C[C@@H](O)[C@@H](C)O"
    result, reason = is_diol(example_smiles)
    print("Is diol?:", result)
    print("Reason:", reason)
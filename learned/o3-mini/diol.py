"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol
A diol is defined here as a compound that contains exactly two alcoholic hydroxy groups.
We ignore –OH groups that are part of carboxylic acids or attached to aromatic carbons.
Examples of diols include (R,R)-butane-2,3-diol, 1,8-tetradecanediol, etc.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol (in our definition) contains exactly 2 alcoholic hydroxy groups.
    Here, an alcoholic hydroxy group is an –OH where the oxygen is bonded to a carbon that:
      – is not aromatic, and
      – is not involved in a carbonyl bond (C=O).
      
    This is done so that carboxylic acid –OH or phenolic –OH groups are not counted.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is (classified as) a diol, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for ANY hydroxyl group (–OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxy_pattern is None:
        return False, "Failed to create hydroxy group pattern"
    
    all_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # We now filter the matches. For each hydroxy group,
    # find the heavy (non-hydrogen) neighbor. We want it to be a carbon,
    # and that carbon must NOT be aromatic and should not be double-bonded to an oxygen.
    count_alcoholic_OH = 0
    for match in all_matches:
        # match is a tuple with the index of the oxygen atom
        o_atom = mol.GetAtomWithIdx(match[0])
        heavy_neighbors = [a for a in o_atom.GetNeighbors() if a.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            # unusual – skip if not exactly one heavy neighbor
            continue
        neigh = heavy_neighbors[0]
        # Only count if neighbor is carbon.
        if neigh.GetAtomicNum() != 6:
            continue
        # Exclude if the carbon is aromatic (i.e. phenol)
        if neigh.GetIsAromatic():
            continue
        # Exclude if the carbon is (likely) part of a carbonyl.
        has_carbonyl = False
        for bond in neigh.GetBonds():
            # Look at bonds from the heavy neighbor; if any bond is a double bond to an oxygen, mark as carbonyl.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(neigh)
                if other.GetAtomicNum() == 8:
                    has_carbonyl = True
                    break
        if has_carbonyl:
            continue
        
        # Otherwise, count this –OH as an alcoholic hydroxyl
        count_alcoholic_OH += 1

    # Now classify as diol only if exactly two alcoholic –OH groups were found.
    if count_alcoholic_OH == 2:
        return True, "Contains exactly two alcoholic hydroxy groups (diol)"
    else:
        return False, f"Found {count_alcoholic_OH} alcoholic hydroxy group(s); a diol requires exactly 2"

# Example usage: (you can remove or update this section if only the function is needed)
if __name__ == "__main__":
    # Example diol: (R,R)-butane-2,3-diol
    example_smiles = "C[C@@H](O)[C@@H](C)O"
    result, reason = is_diol(example_smiles)
    print("Is diol?:", result)
    print("Reason:", reason)
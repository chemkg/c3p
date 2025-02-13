"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
Definition: A fatty acid with a terminal carboxylic acid group (-C(=O)O)
and a hydroxyl (-OH) group attached at the beta (3-) carbon (i.e. two bonds away from the acid carbon).
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid must have a terminal carboxylic acid group and, relative to a chosen acid group,
    an –OH group attached on the beta carbon (the carbon two bonds away from the acid carbon).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens to ensure hydroxyl groups are represented.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for a carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Check for a beta carbon (two bonds away from acid carbon) with an attached hydroxyl group.
    # We use a shortest path approach: the acid carbon (path length 0) to alpha carbon (path length 1)
    # and then beta carbon (path length 2). For each carboxylic acid match, we try to identify such a beta carbon.
    for match in acid_matches:
        acid_idx = match[0]  # the acid carbon from the acid group SMARTS
        # For each carbon atom (atomic num 6) in the molecule except the acid carbon:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetIdx() == acid_idx:
                continue
            # Calculate the shortest path from the acid carbon to this carbon.
            path = Chem.GetShortestPath(mol, acid_idx, atom.GetIdx())
            # A valid beta position should be exactly two bonds away (i.e. 3 atoms long: acid -> alpha -> beta)
            if len(path) != 3:
                continue
            beta_atom = atom  # candidate beta carbon
            # Check if beta_atom carries a hydroxyl group: one of its neighbors must be oxygen with at least one hydrogen.
            for nbr in beta_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check bond type is single bond
                    bond = mol.GetBondBetweenAtoms(beta_atom.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue
                    # If the oxygen atom has at least one hydrogen attached, we consider it a hydroxyl.
                    if nbr.GetTotalNumHs() >= 1:
                        return True, "Found beta-hydroxy (3-) group attached to the fatty acid chain"
    # If no acid group yielded a beta carbon with a hydroxyl group, classification fails.
    return False, "No beta-hydroxy group found on the fatty acid chain"

# Example testing (uncomment the lines below to run a simple test)
# test_smiles = "OCCCCCCCCCCCC[C@@H](O)CC(O)=O"  # (3R)-3,15-dihydroxypentadecanoic acid
# result, reason = is_3_hydroxy_fatty_acid(test_smiles)
# print(result, reason)
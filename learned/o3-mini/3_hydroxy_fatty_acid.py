"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
Definition: Any fatty acid with a terminal carboxylic acid group (-C(=O)O) and a hydroxy (-OH)
functional group attached at the beta (3-) carbon (i.e. the carbon two bonds away from the acid carbon).
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    To be classified as such, the molecule must have a terminal carboxylic acid group and an –OH
    group attached to the beta (3-) carbon of the fatty acid chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure -OH groups are properly represented.
    mol = Chem.AddHs(mol)

    # Define SMARTS for a carboxylic acid group.
    # This looks for a carbon with a doubly bonded oxygen and a single-bonded hydroxyl oxygen.
    acid_smarts = "[CX3](=O)[OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Check each acid group match to see if the beta carbon (second carbon from the acid carbon)
    # has an attached hydroxyl (-OH) group.
    for match in acid_matches:
        # In the SMARTS, match[0] is the carboxylic acid carbon.
        acid_carbon = mol.GetAtomWithIdx(match[0])
        
        # Identify the alpha carbon: a carbon directly connected to the acid carbon.
        alpha_candidates = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not alpha_candidates:
            continue  # Try another acid group if no alpha carbon is found.
        # In a linear fatty acid we assume the first found alpha carbon is correct.
        alpha_atom = alpha_candidates[0]
        
        # The beta carbon is a neighbor of the alpha carbon that is not the acid carbon.
        beta_candidates = [nbr for nbr in alpha_atom.GetNeighbors() 
                           if nbr.GetIdx() != acid_carbon.GetIdx() and nbr.GetAtomicNum() == 6]
        if not beta_candidates:
            continue  # If no beta candidate, move on.
        beta_atom = beta_candidates[0]
        
        # Check if the beta carbon has a hydroxyl (-OH) group attached.
        # That means it must be bonded to an oxygen with at least one hydrogen (explicitly present).
        for nbr in beta_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(beta_atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if nbr.GetTotalNumHs() >= 1:
                        return True, "Found an –OH group on the beta (3-) carbon relative to the acid group"
    
    # If we exit the loop without finding a suitable beta hydroxyl.
    return False, "No beta-hydroxy group found on the fatty acid chain"
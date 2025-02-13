"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
Definition: A fatty acid with a terminal carboxylic acid group (-C(=O)O)
and a hydroxyl (-OH) group attached at the beta (3-) carbon (i.e. the carbon two bonds away from the acid carbon).
This implementation uses an explicit stepwise chain detection starting from the acid group.
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    This implementation first confirms the presence of a terminal carboxylic acid.
    Then, it finds (assuming a linear fatty acid chain) the alpha carbon (neighbor to acid carbon)
    and then the beta (3-) carbon (neighbor of alpha not being the acid carbon).
    Finally, it checks that this beta carbon is substituted with a hydroxyl group (-OH).
    
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
    
    # Add explicit hydrogens to ensure hydroxyl groups are represented correctly.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for a carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H1]"  # Matches a carboxyl group
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Loop over acid matches; often there should be one terminal acid group.
    for match in acid_matches:
        acid_idx = match[0]  # acid carbon (the C of the -C(=O)OH)
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        
        # Identify the alpha carbon: the acid carbon should be bonded to one carbon that is not oxygen.
        alpha_candidates = []
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                alpha_candidates.append(nbr)
        if len(alpha_candidates) != 1:
            # Either chain is not linear or ambiguous; move to next acid match if possible.
            continue
        alpha_atom = alpha_candidates[0]
        
        # Identify the beta (3-) carbon: from alpha, there must be one carbon neighbor that is not the acid carbon.
        beta_candidates = []
        for nbr in alpha_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() == acid_idx:  # skip backwards
                continue
            beta_candidates.append(nbr)
        if len(beta_candidates) != 1:
            # If there is branching or no clear beta, skip this acid center.
            continue
        beta_atom = beta_candidates[0]
        
        # Check that the beta carbon carries a hydroxyl group.
        # We look for an oxygen neighbor bound via a single bond and having at least one explicit hydrogen.
        hydroxyl_found = False
        for nbr in beta_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(beta_atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # Check if the oxygen has at least one hydrogen attached.
                # Note: using GetTotalNumHs() counts explicit and implicit hydrogens.
                if nbr.GetTotalNumHs() >= 1:
                    hydroxyl_found = True
                    break
        
        if hydroxyl_found:
            return True, "Found beta-hydroxy (3-) group attached to the fatty acid chain"
    
    # If none of the acid groups yielded a valid beta-hydroxy, report failure.
    return False, "No beta-hydroxy group found on the fatty acid chain"

# Example testing (uncomment the lines below to run a simple test)
# test_smiles = "OCCCCCCCCCCCC[C@@H](O)CC(O)=O"  # (3R)-3,15-dihydroxypentadecanoic acid
# result, reason = is_3_hydroxy_fatty_acid(test_smiles)
# print(result, reason)
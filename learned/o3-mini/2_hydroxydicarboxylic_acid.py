"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid

Definition: A dicarboxylic acid carrying a hydroxy group on the carbon atom at the position alpha
to (at least one) carboxy group.
"""

from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on the given SMILES.
    The criteria are:
        1. The molecule must contain at least 2 carboxy acid groups (-C(=O)OH).
        2. At least one of the carboxy acid groups must have an alpha-carbon (the carbon attached to
           the carboxyl carbon that is not part of the -C(=O)OH itself) which carries a hydroxy (-OH) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule satisfies the criteria, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string to create a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that -OH groups become explicit
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a carboxy acid group: -C(=O)OH 
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Count unique carboxyl carbon atoms (first atom in our SMARTS match)
    acid_carbons = set(match[0] for match in acid_matches)
    if len(acid_carbons) < 2:
        return False, "Less than 2 carboxy acid groups found"
    
    # Now, for each carboxyl acid group, find the neighbor (the alpha-carbon) that is not
    # an oxygen (part of the carboxyl group) and see if it bears an -OH group.
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        # Get all neighbors of this acid carbon
        for neighbor in acid_carbon.GetNeighbors():
            # We skip oxygen atoms; the alpha-carbon is usually carbon
            if neighbor.GetAtomicNum() != 6:
                continue
            # For this potential alpha-carbon, check if it has an -OH substituent.
            # (Exclude the acid carbon itself when checking for -OH groups.)
            for sub_neigh in neighbor.GetNeighbors():
                if sub_neigh.GetIdx() == acid_carbon_idx:
                    continue
                # Check if the neighbor is an oxygen atom
                if sub_neigh.GetAtomicNum() == 8:
                    # Check if this oxygen is bonded to at least one hydrogen (i.e. an -OH)
                    has_hydrogen = any(n.GetAtomicNum() == 1 for n in sub_neigh.GetNeighbors())
                    if has_hydrogen:
                        return True, "Found a carboxy group with an alpha-carbon bearing a hydroxy group"
                        
    return False, "No carboxyl group found with an alpha-carbon that has a hydroxy substituent"

# For testing, you can uncomment lines below and run various examples:
# examples = [
#     "OC(=O)C(CCC(O)=O)O",  # 2-hydroxyadipic acid example (should match)
#     "O=C(O)CC(O)=O",       # Does not have an alpha-hydroxy (example dicarboxylic acid)
# ]
# for smi in examples:
#     result, reason = is_2_hydroxydicarboxylic_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")
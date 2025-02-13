"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated fatty acid.
Definition: Any fatty acid with one (and only one) carbon–carbon unsaturation 
(double or triple bond) in the acyl chain and a free (terminal) carboxylic acid group.
MUFAs are known to have positive effects on the cardiovascular system and in diabetes treatment.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    
    The criteria applied here are:
      1) The molecule must have exactly one carboxylic acid group. We use a SMARTS pattern
         that covers both protonated (-C(=O)[OH]) and deprotonated (-C(=O)[O-]) forms.
      2) The acid group must be terminal (its carboxyl carbon is attached to exactly one carbon atom).
      3) The carbon chain attached to the acid group ("acyl chain") is extracted (only traversing carbon atoms).
      4) Within that acyl chain, there must be exactly one unsaturation (a carbon–carbon double or triple bond).
      5) We also require that the acyl chain is not trivially short (at least 3 carbon atoms beyond the acid, 
         though this threshold can be adjusted).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets all criteria for being a monounsaturated fatty acid, 
              False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a carboxylic acid group: match protonated or deprotonated forms.
    # [CX3](=O)[OX2H1] matches a protonated acid, and [CX3](=O)[OX1H0] will match 
    # a negatively charged carboxylate oxygen. We combine them into a single pattern:
    acid_smarts = "[CX3](=O)[OX1H0,OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        # Either no acid group or multiple acid groups were found.
        return False, f"Found {len(acid_matches)} acid group(s); need exactly one terminal carboxylic acid group"
        
    # For the (only) matching acid group, identify the carboxyl carbon.
    # By our SMARTS, the first atom is the carboxyl carbon.
    acid_match = acid_matches[0]
    acid_carbon_idx = acid_match[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Check that the acid carbon is "terminal": it should be bonded to exactly one carbon atom.
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Acid group is not terminal (its acid carbon is attached to more than one carbon)"
    
    # The only carbon neighbor is the alpha (first) carbon of the fatty acyl chain.
    alpha = carbon_neighbors[0]
    
    # Now, extract the contiguous acyl chain: all carbon atoms connected to the alpha carbon 
    # while not traversing back into the acid group.
    chain_atoms = set()
    def dfs(atom, coming_from):
        if atom.GetIdx() in chain_atoms:
            return
        chain_atoms.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            # Only traverse if the neighbor is a carbon and is not the atom we came from.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != coming_from:
                dfs(nbr, atom.GetIdx())
    dfs(alpha, acid_carbon_idx)
    
    # For a fatty acid, the acyl chain should be non-trivial.
    if len(chain_atoms) < 3:
        return False, f"Acyl chain too short (only {len(chain_atoms)} carbon(s) found beyond the acid group)"
    
    # Now count carbon–carbon unsaturation bonds (double or triple) within the acyl chain.
    unsaturation_count = 0
    # For every bond in the molecule, if both atoms are in the chain_atoms set then
    # check if that bond is a double or triple bond.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        if a1 in chain_atoms and a2 in chain_atoms:
            btype = bond.GetBondType()
            if btype in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
                unsaturation_count += 1
                
    if unsaturation_count != 1:
        return False, f"Found {unsaturation_count} carbon–carbon unsaturation(s) in the acyl chain; need exactly one"
    
    return True, "Molecule contains one terminal carboxylic acid group and an acyl chain with exactly one carbon–carbon unsaturation, qualifying as a monounsaturated fatty acid"

# Uncomment the lines below to run a quick test on one of the provided examples.
# example_smiles = "CCCCCC\\C=C\\CCCCCCCCCC(O)=O"  # trans-vaccenic acid
# result, reason = is_monounsaturated_fatty_acid(example_smiles)
# print(result, reason)
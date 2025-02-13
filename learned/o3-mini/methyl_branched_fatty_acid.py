"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
Definition: any branched-chain fatty acid containing methyl branches only.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a methyl-branched fatty acid.
    A methyl-branched fatty acid must:
      - Contain a carboxylic acid group (COOH or its deprotonated form).
      - Have a long alkyl chain (starting from the alpha carbon attached to the carboxyl group).
      - Possess branches on the alkyl chain and every branch must be an isolated methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES and add hydrogens (so that we can check H count on atoms)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for a carboxylic acid group (covers COOH and its deprotonated form)
    acid_smarts = "[CX3](=O)[OX2H1,OX1-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # For our purposes, take the first carboxyl match.
    # In the SMARTS, the first atom should be the carboxyl carbon.
    carboxyl_idx = acid_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Find the alpha carbon: a neighbor (atomic number 6) of the carboxyl carbon.
    alpha_atom = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # carbon
            alpha_atom = nbr
            break
    if alpha_atom is None:
        return False, "No alpha carbon found attached to carboxyl group"
    
    # Now, obtain the longest carbon chain starting from the alpha carbon.
    # We perform a depth-first search (DFS) on carbon atoms only.
    def dfs(atom, visited):
        visited = visited + [atom.GetIdx()]
        # Find all carbon neighbors not yet visited.
        paths = []
        extendable = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited:
                continue
            extendable = True
            new_path = dfs(nbr, visited)
            paths.append(new_path)
        if not extendable:
            return visited
        else:
            # Return the longest path (by number of atoms)
            best = max(paths, key=lambda p: len(p))
            return best

    main_chain = dfs(alpha_atom, [])
    
    if len(main_chain) < 2:
        return False, "The alkyl chain is too short"
    
    # Now, examine branches: for every atom in the main chain examing neighbors that are carbon and not in the main chain.
    branch_count = 0
    for idx in main_chain:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in main_chain:
                continue
            # Check if this neighbor is a methyl group:
            # A methyl group should have exactly one heavy-atom neighbor (the main-chain carbon)
            heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1]
            # Using the hydrogens we added, the methyl carbon should have 1 heavy neighbor and 3 hydrogens.
            if len(heavy_neighbors) != 1 or nbr.GetTotalNumHs() != 3:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that is not a simple methyl group"
            branch_count += 1
    
    if branch_count == 0:
        return False, "No methyl branches detected on the fatty acid chain"
    
    return True, "Molecule contains a carboxyl group with an alkyl chain having only methyl branches"

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "OC(=O)CC(C)=C",  # Isopropenylacetic acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methyltriacontanoic acid
        "CCCCCCCCCCCCCCCCCCCCCC(O)=O"  # a linear fatty acid (should fail)
    ]
    for smi in test_smiles:
        result, reason = is_methyl_branched_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
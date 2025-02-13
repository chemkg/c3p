"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
Definition: any branched-chain fatty acid containing methyl branches only.
The molecule must have a carboxylic acid group and an alkyl chain (starting at the alpha-carbon)
whose only branches are methyl groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a methyl-branched fatty acid.
    
    Requirements:
      - Contains a carboxylic acid group (protonated or deprotonated).
      - Has an alkyl chain starting from the alpha-carbon that is longest in the carbon-only subgraph.
      - Any branch (carbon substituent off the main chain) MUST be a simple methyl group:
          That means the substituent carbon has only one heavy-atom neighbor (the main-chain carbon)
          and carries exactly 3 hydrogens.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES and add explicit hydrogens for proper H count checks.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)

    # Define SMARTS for a carboxylic acid (covers COOH and its deprotonated form)
    acid_smarts = "[CX3](=O)[OX2H1,OX1-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # For our purposes, take the first carboxyl match.
    # The first atom in the match is the carboxyl carbon.
    carboxyl_idx = acid_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Find the alpha carbon: a neighbor carbon (atomic number 6) of the carboxyl carbon.
    alpha_atom = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            alpha_atom = nbr
            break
    if alpha_atom is None:
        return False, "No alpha carbon found attached to carboxyl group"
    
    # Build the carbon-only subgraph starting from the alpha carbon.
    # We want to obtain the longest continuous chain (as a list of atom indices).
    # Exclude the carboxyl carbon from being included, so that the chain
    # does not “backtrack” into the acid head.
    def dfs(current_atom, visited):
        visited = visited + [current_atom.GetIdx()]
        extended_paths = []
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue  # only consider carbon atoms
            # Do not go into forbidden atoms (carboxyl carbon) or already visited atoms.
            if nbr.GetIdx() in visited or nbr.GetIdx() == carboxyl_idx:
                continue
            new_path = dfs(nbr, visited)
            extended_paths.append(new_path)
        if not extended_paths:
            return visited
        else:
            # Return the longest path found.
            return max(extended_paths, key=lambda p: len(p))
    
    main_chain = dfs(alpha_atom, [])
    if len(main_chain) < 2:
        return False, "The alkyl chain (from the alpha carbon) is too short"
    
    # Now, for every carbon atom in the main chain, examine its carbon neighbors that are NOT in the main chain.
    branch_count = 0
    for idx in main_chain:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in main_chain:
                continue  # neighbor is part of main chain
            # Check if the neighbor is a simple methyl group:
            # A proper methyl group must have exactly one heavy neighbor (the attaching main-chain atom)
            # and exactly 3 bound hydrogens.
            heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1]
            if len(heavy_neighbors) != 1:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that is not a simple methyl group (unexpected heavy neighbor count)"
            # Check total number of hydrogens on this branch carbon.
            if nbr.GetTotalNumHs() != 3:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that does not have exactly 3 hydrogens"
            branch_count += 1

    # If no branches are found, then it is a linear fatty acid (not methyl-branched).
    if branch_count == 0:
        return False, "No methyl branches detected on the fatty acid chain"

    return True, "Molecule has a carboxyl group with an alkyl chain that only has methyl branches"

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "OC(=O)CC(C)=C",  # Isopropenylacetic acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methyltriacontanoic acid
        "CCC(C)C(O)=O",  # 2-methylbutyric acid
    ]
    for smi in test_smiles:
        result, reason = is_methyl_branched_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
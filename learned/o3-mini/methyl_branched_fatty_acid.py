"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
Definition: any branched-chain fatty acid containing methyl branches only.
The molecule must have a carboxylic acid group and an alkyl chain (starting at the alpha-carbon)
whose only branches (i.e. substituents off the main chain) are simple methyl groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule (given by a SMILES string) is a methyl-branched fatty acid.
    
    Requirements:
      - Contains a carboxylic acid (or its deprotonated form) group.
      - Has an alkyl chain starting from the alpha-carbon (attached to the carboxyl carbon).
      - Any branch (carbon substituent off the main chain) MUST be a simple methyl group:
          That means the branch carbon must be sp³, have exactly one heavy-atom neighbor (the main-chain atom)
          and exactly 3 hydrogens.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a methyl-branched fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES; if invalid, return error.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens to properly count them.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for a carboxylic acid group (covers COOH and the deprotonated form).
    acid_smarts = "[CX3](=O)[OX2H1,OX1-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Take the first carboxyl match; the first atom in the match is the carboxyl carbon.
    carboxyl_idx = acid_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Find the alpha-carbon: a neighbor carbon (atomic number 6) of the carboxyl carbon.
    alpha_atom = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            alpha_atom = nbr
            break
    if alpha_atom is None:
        return False, "No alpha carbon found attached to carboxyl group"
    
    # Build the main carbon chain starting from the alpha-carbon.
    # We will perform a DFS on the graph of carbon atoms (atomic num 6) – excluding the carboxyl carbon,
    # to pick the longest chain.
    def dfs(current_atom, visited):
        current_idx = current_atom.GetIdx()
        visited = visited + [current_idx]
        possible_paths = []
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue  # only consider carbons
            if nbr.GetIdx() == carboxyl_idx:
                continue  # do not go back into the acid head
            if nbr.GetIdx() in visited:
                continue
            path = dfs(nbr, visited)
            possible_paths.append(path)
        if not possible_paths:
            return visited
        # return the longest path found
        return max(possible_paths, key=lambda p: len(p))
    
    main_chain = dfs(alpha_atom, [])
    if len(main_chain) < 2:
        return False, "The alkyl chain (starting at the alpha carbon) is too short"
    
    # For every atom in the main chain, check for branches that are not part of the main chain.
    branch_count = 0
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only consider carbon neighbors not in the main chain.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in main_chain:
                continue  # neighbor is part of main chain
            # Now, check that this branch substituent is a simple methyl group.
            # a simple methyl should be sp3 hybridized.
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that is not sp3 hybridized (likely not a simple methyl group)"
            # Count heavy (non-H) neighbors of the branch atom. It should have exactly 1 (the attachment).
            heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1]
            if len(heavy_neighbors) != 1:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that is not a simple methyl group (unexpected heavy neighbor count)"
            # Check that the branch carbon carries exactly 3 hydrogens.
            if nbr.GetTotalNumHs() != 3:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that does not have exactly 3 hydrogens"
            branch_count += 1

    # If no branches were found, it is a linear fatty acid.
    if branch_count == 0:
        return False, "No methyl branches detected on the fatty acid chain"
    
    return True, "Molecule has a carboxyl group with an alkyl chain that only has methyl branches"

# Example usage:
if __name__ == '__main__':
    test_smiles_list = [
        "OC(=O)CC(C)=C",  # Isopropenylacetic acid (should be methyl-branched)
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methyltriacontanoic acid
        "CCC(C)CCCCCCCCCCCCC(O)=O",  # 14-methylhexadecanoic acid
        "CC(C)CCCCCCCCCC(O)=O",  # isotridecanoic acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 24-methylpentacosanoic acid
        "OC(=O)C(CCC)CC",  # alpha-ethyl valeric acid (should fail since branch is ethyl)
        "CCC(C)CCCCCCCCCCCCCCCCCCC(O)=O",  # 20-methyldocosanoic acid
        "[O-]C(=O)CCC([N+](C)(C)C)C",  # 4-aminovaleric acid betaine
        "CCC(C)C(O)=O",  # 2-methylbutyric acid
    ]
    
    for smi in test_smiles_list:
        result, reason = is_methyl_branched_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: short-chain fatty acid
Definition:
  An aliphatic monocarboxylic acid with a chain length (including the carboxyl carbon) of less than C6.
  In addition, substituents on the backbone must be hydrocarbon in nature – if any non‐hydrocarbon substituent is present
  (apart from the carboxyl group itself), the compound is not normally regarded as a short‐chain fatty acid.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    
    A short-chain fatty acid in this context is defined as an aliphatic monocarboxylic acid
    (i.e. containing exactly one carboxyl group) whose longest continuous carbon chain (including
    the carboxyl carbon) has fewer than 6 carbon atoms. Any substituent on the chain that is not a 
    hydrocarbon (i.e. not a carbon or an implicit hydrogen) would disqualify it.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is acyclic (no rings) because short‐chain fatty acids are open‐chain.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not purely aliphatic"
    
    # Look for a carboxyl group.
    # We use a SMARTS pattern for a carboxylic acid: [C](=O)[O;H1]
    acid_smarts = "[C](=O)[O;H1]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    matches = mol.GetSubstructMatches(acid_query)
    if not matches:
        return False, "No carboxyl (acid) group found"
    if len(matches) != 1:
        return False, "More than one carboxyl group was found, not a monocarboxylic acid"
    
    # The carboxyl carbon is the first atom in the match.
    acid_idx = matches[0][0]
    acid_atom = mol.GetAtomWithIdx(acid_idx)
    if acid_atom.GetAtomicNum() != 6:
        return False, "Carboxyl group does not have a carbon as expected"
    
    # Define a helper function to build the connected "acyl chain" graph.
    # The acyl chain includes carbons connected via C–C bonds from the acid carbon.
    def get_carbon_chain(start_idx):
        chain = set()
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in chain:
                continue
            chain.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # only consider carbon neighbors
                    # Avoid traveling into the acid group oxygens (or other heteroatoms)
                    stack.append(neighbor.GetIdx())
        return chain

    # Build the carbon chain starting at the acid carbon.
    chain_indices = get_carbon_chain(acid_idx)
    
    # Note: the carboxyl carbon itself is part of a functional group,
    # so when looking for non-hydrocarbon substituents, we will ignore its attached oxygen atoms.
    #
    # Now we build an adjacency map for the acyl chain (only carbons within chain_indices)
    adjacency = {idx: [] for idx in chain_indices}
    for idx in chain_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            nidx = neighbor.GetIdx()
            if nidx in chain_indices:
                adjacency[idx].append(nidx)
    
    # Compute the longest simple path in the carbon chain graph.
    # Because short-chain fatty acids have few carbons, a brute-force DFS is acceptable.
    max_path_length = 0  # maximum number of carbon atoms in a simple path

    def dfs(current, visited):
        nonlocal max_path_length
        visited.add(current)
        # update max_path_length based on visited set size
        if len(visited) > max_path_length:
            max_path_length = len(visited)
        for neigh in adjacency[current]:
            if neigh not in visited:
                dfs(neigh, visited)
        visited.remove(current)
    
    dfs(acid_idx, set())

    # Check chain length: by definition the total number of carbons (including the carboxyl carbon)
    # should be less than 6.
    if max_path_length >= 6:
        return False, f"Longest carbon chain has {max_path_length} carbons, exceeding the short-chain limit (<6)"
    
    # Now verify that any substituents on the acyl chain (attached to carbons in chain_indices)
    # that are not part of the acyl chain or the carboxyl acid group are hydrocarbon in nature.
    # Allowed atoms: carbon (atomic number 6) and (implicit) hydrogens.
    # For the acid carbon itself, ignore the two oxygens that form the carboxyl group.
    for idx in chain_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            nidx = neighbor.GetIdx()
            # If neighbor is already part of the carbon chain, skip.
            if nidx in chain_indices:
                continue
            # For the acid carbon, allow the oxygen atoms that form the acid group.
            if idx == acid_idx and neighbor.GetAtomicNum() == 8:
                continue
            # If any neighbor is not a hydrogen (atomic number 1) or carbon, disqualify.
            if neighbor.GetAtomicNum() not in (1, 6):
                return False, "Non-hydrocarbon substituent found attached to the acyl chain"
    
    return True, "Molecule is an aliphatic monocarboxylic acid with a carbon chain length of less than C6"

# (Optional) You can include a main block to test several SMILES strings.
if __name__ == "__main__":
    test_smiles = [
        "CCCCC(O)=O",         # valeric acid (C5 acid)
        "CC(C)=CC(O)=O",       # 3-methylbut-2-enoic acid (backbone length = 5)
        "CCC(C)C(O)=O",        # (R)-2-methylbutyric acid (backbone length = 4)
        "CCCC(O)=O",           # butyric acid (C4 acid)
        "[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C"  # heliosupine (complex cyclic)
    ]
    
    for smi in test_smiles:
        result, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\n  Is short-chain fatty acid? {result}\n  Reason: {reason}\n")
"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: Octadecadienoic acid
Definition: Any straight-chain, C18 polyunsaturated fatty acid having exactly two C=C bonds.
This improved implementation:
  - Searches for a carboxylic acid group (C(=O)[OH] or O=C(O))
  - Identifies the acid carbon (the carbonyl carbon)
  - Searches via DFS for an unbranched sequence of 18 carbons (including the acid carbon) 
    connected only by single or double bonds (none of which are in rings)
  - Verifies that in that chain exactly 2 bonds are double bonds.
Note: This is a simplified approach that tries to prevent false positives by rejecting chains that have extra carbon attachments.
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines whether a molecule is an octadecadienoic acid based on its SMILES string.
    The criteria are:
      - A carboxylic acid group is present.
      - One can identify a straight, unbranched sequence (path) of 18 carbon atoms (including the acid carbon)
        from the acid carbon to a terminal methyl group.
      - All bonds along this path are single or double and none is part of a ring.
      - Exactly 2 of these bonds are double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as an octadecadienoic acid, False otherwise.
        str: A textual explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Look for a carboxylic acid group.
    # We search for the common representations of acid: either "C(=O)[OH]" or "O=C(O)"
    acid_smarts_list = ["C(=O)[OH]", "O=C(O)"]
    acid_matches = []
    acid_pattern = None
    for smarts in acid_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            acid_matches = matches
            acid_pattern = pattern
            break
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Assume the first match is our acid group.
    # Identify the acid carbon (the carbon with atomic number 6 among the matched atoms)
    acid_carbon = None
    for atom_idx in acid_matches[0]:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            acid_carbon = atom_idx
            break
    if acid_carbon is None:
        return False, "Acid group found but no acid carbon identified"
    
    # Precompute candidate terminal methyl carbons.
    # A terminal methyl (at the other end of the chain) should have exactly one carbon neighbor.
    candidate_terminals = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() != acid_carbon:
            carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 1:
                candidate_terminals.append(atom.GetIdx())
    if not candidate_terminals:
        return False, "No candidate terminal methyl group found"
    
    # We want a path of 18 carbon atoms (chain length) from acid_carbon to one terminal.
    target_chain_len = 18  # number of carbons
    valid_chain_found = False
    found_details = []
    
    # Define DFS that only follows bonds that are SINGLE or DOUBLE and are not in rings.
    def dfs_paths(current, target, max_len, path, paths):
        if len(path) == max_len:
            if current == target:
                paths.append(path.copy())
            return
        for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
            # Only follow carbon atoms.
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in path:
                continue
            bond = mol.GetBondBetweenAtoms(current, nbr_idx)
            # Exclude bonds that are not single or double or if they are part of a ring.
            if bond is None:
                continue
            if bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                continue
            if bond.IsInRing():
                continue
            path.append(nbr_idx)
            dfs_paths(nbr_idx, target, max_len, path, paths)
            path.pop()
    
    # Try every candidate terminal methyl group.
    for term in candidate_terminals:
        paths = []
        dfs_paths(acid_carbon, term, target_chain_len, [acid_carbon], paths)
        for p in paths:
            # Verify path has exactly target_chain_len carbons.
            if len(p) != target_chain_len:
                continue

            # Check that the chain is linear:
            # a truly "straight" chain should have no additional carbon (branch) attachments:
            branch_found = False
            for idx in p:
                atom = mol.GetAtomWithIdx(idx)
                # among carbon neighbors in the molecule, all should be part of the chain.
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in p:
                        branch_found = True
                        break
                if branch_found:
                    break
            if branch_found:
                found_details.append(f"Chain {p} is branched (extra carbon neighbor found).")
                continue

            # Count double bonds along consecutive atoms in the chain.
            dbl_bond_count = 0
            valid_bonds = True
            for i in range(len(p)-1):
                bond = mol.GetBondBetweenAtoms(p[i], p[i+1])
                if bond is None:
                    valid_bonds = False
                    break
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    dbl_bond_count += 1
                # (Bond type was already restricted in DFS so no further check here)
            if not valid_bonds:
                found_details.append(f"Chain {p} has an invalid bond between consecutive atoms.")
                continue
            
            if dbl_bond_count != 2:
                found_details.append(f"Chain {p} has {dbl_bond_count} C=C bond(s), expected exactly 2.")
                continue

            # If a valid chain is found, we are done.
            valid_chain_found = True
            break
        if valid_chain_found:
            break
    
    if valid_chain_found:
        return True, "Molecule is a straight-chain C18 fatty acid with exactly 2 C=C bonds"
    else:
        debug_reasons = "; ".join(found_details) if found_details else "No valid chain with required criteria found."
        return False, f"No valid straight-chain C18 found with exactly 2 C=C bonds. {debug_reasons}"

# Example usage:
if __name__ == "__main__":
    # Testing one provided example: 9-cis,11-trans-octadecadienoic acid
    test_smiles = "CCCCCC\\C=C\\C=C/CCCCCCCC(O)=O"
    result, reason = is_octadecadienoic_acid(test_smiles)
    print(result, reason)
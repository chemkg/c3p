"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: Octadecadienoic acid
Definition: Any straight-chain, C18 polyunsaturated fatty acid having exactly two C=C double bonds.
This implementation:
  - Searches for a carboxylic acid group (C(=O)[OH])
  - Selects the acid carbon and then attempts to find a simple (acyclic) carbon-only path
    from that atom to a terminal methyl group that contains exactly 18 carbons.
  - In that candidate chain, every bond must be either a single or double bond (no triple bonds)
  - And exactly two of the chain bonds must be carbon–carbon double bonds.
Note: This is a simplified approach and may not treat all edge cases.
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is defined as any straight-chain, C18 fatty acid that has exactly
    two C=C double bonds along the chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies, False otherwise.
        str: A textual reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Look for the carboxylic acid group using a SMARTS pattern.
    # We require C(=O)[OH] where the acid carbon is the one with atomic number 6.
    acid_smarts = "C(=O)[OH]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Assume the first match is our acid group.
    acid_carbon = None
    for atom_idx in acid_matches[0]:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            acid_carbon = atom_idx
            break
    if acid_carbon is None:
        return False, "Carboxylic acid group found but no acid carbon identified"
    
    # Helper: DFS to recursively find carbon-only paths of a given length (number of atoms).
    # We are only following carbon atoms (atomic number 6).
    def dfs_paths(current, target, max_len, path, paths):
        if len(path) == max_len:
            if current == target:
                paths.append(path.copy())
            return
        # Traverse neighbors that are carbons and not already in the path.
        for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in path:
                continue
            path.append(nbr.GetIdx())
            dfs_paths(nbr.GetIdx(), target, max_len, path, paths)
            path.pop()
    
    # Identify candidate terminal methyl carbons.
    # For our purposes a terminal methyl carbon must be carbon with exactly 1 carbon neighbor.
    candidate_terminals = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Count neighbors that are carbon atoms.
        carbon_neighs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighs) == 1 and atom.GetIdx() != acid_carbon:
            candidate_terminals.append(atom.GetIdx())
    if not candidate_terminals:
        return False, "No terminal methyl group found (as a candidate end of chain)"
    
    # We now search for a carbon-only path from the acid carbon to a terminal carbon.
    # In a C18 fatty acid the chain should have exactly 18 carbon atoms (including the acid carbon)
    target_chain_len = 18
    valid_chain_found = False
    found_details = []
    
    for term in candidate_terminals:
        paths = []
        dfs_paths(acid_carbon, term, target_chain_len, [acid_carbon], paths)
        for p in paths:
            # p is a candidate chain: it must be exactly 18 carbons.
            if len(p) != target_chain_len:
                continue
            # Check that the path is “straight”.
            # For each internal carbon (other than first and last), the number of adjacent carbons that are in the chain should be exactly 2.
            straight = True
            for i, idx in enumerate(p):
                # Get neighbors in chain
                chain_neighs = [nbr for nbr in mol.GetAtomWithIdx(idx).GetNeighbors() 
                                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in p]
                if i == 0 or i == len(p)-1:
                    if len(chain_neighs) != 1:
                        straight = False
                        break
                else:
                    if len(chain_neighs) != 2:
                        straight = False
                        break
            if not straight:
                found_details.append(f"Chain {p} is branched.")
                continue

            # Check the bonds between consecutive carbons in the chain:
            dbl_bond_count = 0
            valid_bond_types = True
            for i in range(len(p)-1):
                bond = mol.GetBondBetweenAtoms(p[i], p[i+1])
                if bond is None:
                    valid_bond_types = False
                    break
                btype = bond.GetBondType()
                # Allow only SINGLE or DOUBLE bonds.
                if btype not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                    valid_bond_types = False
                    break
                if btype == Chem.BondType.DOUBLE:
                    dbl_bond_count += 1
            if not valid_bond_types:
                found_details.append(f"Chain {p} contains a bond that is not single or double.")
                continue
            
            if dbl_bond_count != 2:
                found_details.append(f"Chain {p} has {dbl_bond_count} C=C bond(s), expected exactly 2.")
                continue

            # If we are here, we found a valid chain.
            valid_chain_found = True
            break
        if valid_chain_found:
            break

    if valid_chain_found:
        return True, "Molecule is a straight-chain C18 fatty acid with exactly 2 C=C bonds"
    else:
        # We could not find any valid path. Report the reasons from candidate chains.
        debug_reasons = "; ".join(found_details) if found_details else "No chain found with the required criteria."
        return False, f"No valid straight-chain C18 found with exactly 2 C=C bonds. {debug_reasons}"

# Example usage:
if __name__ == "__main__":
    # test one of the provided examples:
    test_smiles = "CCCCCC\\C=C\\C=C/CCCCCCCC(O)=O"  # 9-cis,11-trans-octadecadienoic acid
    result, reason = is_octadecadienoic_acid(test_smiles)
    print(result, reason)
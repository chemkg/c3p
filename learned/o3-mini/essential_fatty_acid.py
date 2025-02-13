"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential Fatty Acid
Definition: Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
This function attempts a more focused approach:
  1. It checks that the SMILES string can be parsed.
  2. It then looks for either a free carboxylic acid group or a fatty acyl ester group.
     (Many fatty acids occur as acyl chains in esters such as in glycerophospholipids.)
  3. If one is found, it “extracts” a contiguous chain of carbon atoms (ignoring any atoms in rings)
     starting at the carbonyl carbon. Fatty acyl chains are usually linear and unbranched.
  4. It counts the number of carbons in that chain and the number of carbon–carbon double bonds.
  5. If the chain is long enough (heuristically ≥ 12 carbons) and has at least 2 C=C bonds,
     we classify the molecule as an essential fatty acid.
Note: This “chain–extraction” is a heuristic. It builds an induced graph on non‐ring carbons
      and then uses a DFS that does not allow branching (more than one choice) to pick out a linear chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    The approach is:
      1. Verify the molecule is valid.
      2. Look for a free carboxylic acid group or a fatty acyl ester group.
         – For a free acid we use the SMARTS "C(=O)[O;H1,-1]"
         – For an ester we use the SMARTS "C(=O)O[C]" (which will capture fatty acyl units in, e.g., phospholipids)
      3. For each match, take the carboxyl (or acyl) carbon as an “anchor” and try to extract
         a contiguous chain over non‐ring carbons.
      4. Count how many carbon atoms are in that chain and how many double bonds are present.
      5. If the chain is long (≥ 12 carbons) and with at least 2 C=C bonds then we classify as essential fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an essential fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS patterns for a free acid and for an ester fatty acyl carbonyl.
    acid_smarts = "C(=O)[O;H1,-1]"   # carboxyl group (free acid)
    ester_smarts = "C(=O)O[C]"         # acyl ester (e.g., as in phospholipids)
    pat_acid = Chem.MolFromSmarts(acid_smarts)
    pat_ester = Chem.MolFromSmarts(ester_smarts)
    
    # Collect candidate indices (the carboxyl or acyl carbon) from matches.
    candidate_indices = set()
    for match in mol.GetSubstructMatches(pat_acid):
        candidate_indices.add(match[0])
    for match in mol.GetSubstructMatches(pat_ester):
        candidate_indices.add(match[0])
    
    if not candidate_indices:
        return False, "No carboxylic acid (or fatty acyl ester) functional group found."
    
    # Build an induced graph (as a dictionary) of all carbon atoms that are not in any ring.
    # These are the atoms likely to be part of a linear fatty acyl chain.
    allowed = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            idx = atom.GetIdx()
            # For each such carbon, list its neighbors that are also allowed.
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() 
                         if nbr.GetAtomicNum() == 6 and not nbr.IsInRing()]
            allowed[idx] = neighbors
    
    # For each candidate starting atom, attempt to extract a linear (nonbranching) chain.
    best_chain = []
    best_dbonds = 0  # number of C=C bonds in the best chain
    # We use DFS that does not allow a branching (more than one extension from any carbon after the start).
    for start in candidate_indices:
        if start not in allowed:
            continue  # if the anchor is not in the allowed (non‐ring C) set, skip
        paths = []  # will store linear paths (lists of atom indices)
        def dfs(current, parent, path):
            new_path = path + [current]
            # Get allowed carbon neighbors that are not the parent.
            neighs = [n for n in allowed.get(current, []) if n != parent]
            # For a linear chain we expect at most one extension.
            if len(neighs) != 1:
                # End the chain here.
                paths.append(new_path)
            else:
                dfs(neighs[0], current, new_path)
        dfs(start, None, [])
        # Among all paths from this starting point, pick the longest.
        for p in paths:
            if len(p) > len(best_chain):
                best_chain = p
                # Count the number of double bonds along this chain.
                count_db = 0
                for i in range(len(p) - 1):
                    bond = mol.GetBondBetweenAtoms(p[i], p[i+1])
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        count_db += 1
                best_dbonds = count_db

    # If no chain was found, then fail.
    if not best_chain:
        return False, "No linear fatty acyl chain found starting from an acid/ester group."
    
    # Check the chain length and the degree of unsaturation.
    chain_length = len(best_chain)
    if chain_length < 12:
        return False, f"Fatty acyl chain too short ({chain_length} carbons)."
    if best_dbonds < 2:
        return False, f"Insufficient polyunsaturation (found {best_dbonds} C=C bonds, need at least 2)."
    
    return True, (f"Contains fatty acyl chain with {chain_length} carbons "
                  f"and {best_dbonds} C=C bonds, meeting criteria for an essential fatty acid.")

# (Optional) Testing on one of the provided examples:
if __name__ == '__main__':
    test_smiles = [
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCC(O)=O",  # all-cis-7,10,13-hexadecatrienoic acid
        "OC(=O)CC",  # too short example
    ]
    for smi in test_smiles:
        result, reason = is_essential_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
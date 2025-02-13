"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: Short-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation
of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    The algorithm checks:
      - The presence of a thioester group (C(=O)S).
      - That the acyl (fatty acid) fragment attached to the carbonyl carbon is short 
        (we require no more than 6 carbon atoms in that fragment).
      - The presence of a CoA moiety by detecting an adenine fragment using multiple SMARTS.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES input
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a thioester functionality: a carbonyl linked to a sulfur.
    thioester_smarts = "[C](=O)[S]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond detected, hence not an acyl-CoA"
    
    # For this check, consider only the first found thioester bond.
    # The pattern returns a tuple of atom indices: (acyl_carbon, carbonyl_oxygen, sulfur)
    acyl_carbon_idx, _, sulfur_idx = thioester_matches[0]
    
    # Define a recursive depth-first search (DFS) function to count carbon atoms in the acyl fragment.
    # We want to count only carbons reachable from the acyl carbon without crossing the bond to sulfur.
    def dfs_count_atoms(atom_idx, visited):
        count = 0
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            count += 1
        visited.add(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Avoid traversing the bond to the sulfur (connection to CoA part)
            if nbr_idx == sulfur_idx:
                continue
            # traverse only if the neighbor is carbon
            if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                count += dfs_count_atoms(nbr_idx, visited)
        return count

    visited = set()
    acyl_carbons = dfs_count_atoms(acyl_carbon_idx, visited)
    
    # Define the allowed acyl chain length. Here we require a short-chain fatty acid:
    # For instance, acetyl (2 carbons), propionyl (3), butyryl (4), valeryl (5), caproyl (6).
    if acyl_carbons > 6:
        return False, f"Acyl chain contains {acyl_carbons} carbons, too long for a short-chain fatty acid"
    if acyl_carbons < 2:
        return False, "Acyl chain too short (less than 2 carbons) to be a fatty acid"
    
    # Check for a fragment typical of CoA. Many CoA molecules contain an adenine moiety.
    # The exact SMARTS representation of adenine can vary, so we include several variants.
    adenine_smarts_list = [
        "n1cnc2c(n)cnc12",      # original pattern
        "n1c2ncnc2ncn1",        # alternative ordering
        "c1nc2nc[nH]c2n1",      # another common form
    ]
    adenine_found = False
    for smarts in adenine_smarts_list:
        adenine_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(adenine_pattern):
            adenine_found = True
            break

    if not adenine_found:
        return False, "No adenine fragment typical of CoA detected"
    
    # If all checks passed, the molecule satisfies the criteria.
    return True, f"Found thioester bond with an acyl chain of {acyl_carbons} carbons and a CoA moiety"

# Optional testing examples. You can run this module to test one example.
if __name__ == "__main__":
    # Example test: (R)-3-hydroxypentanoyl-CoA from the provided examples.
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C[C@@H](CC)O)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
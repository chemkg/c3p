"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: Short-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the condensation
of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    The algorithm checks:
     - That a thioester (C(=O)S) functional group is present.
     - That the acyl (fatty acid) fragment attached to the carbonyl carbon is 
       "short" (here defined by having no more than 6 carbon atoms in total).
     - That a fragment typical of the CoA moiety, using an adenine ring fragment,
       is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for thioester: acyl (C=O) bound to S.
    thioester_smarts = "[C](=O)[S]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond detected, hence not an acyl-CoA"
    
    # For the purpose of this classification we check the first found thioester.
    # thioester_match is a tuple of atom indices: (acyl_carbon, carbonyl_oxygen, sulfur)
    acyl_carbon_idx, _, sulfur_idx = thioester_matches[0]
    
    # Now, count the number of carbon atoms in the acyl (fatty acid) fragment.
    # We want to count the carbons reachable from the acyl carbon WITHOUT crossing the
    # bond to the sulfur. We perform a DFS that only follows bonds between carbon atoms.
    def dfs_count_atoms(atom_idx, visited):
        count = 0
        atom = mol.GetAtomWithIdx(atom_idx)
        # Count the atom if it is carbon.
        if atom.GetAtomicNum() == 6:
            count += 1
        visited.add(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # do not cross to the sulfur atom (the connection to CoA)
            if nbr_idx == sulfur_idx:
                continue
            # Only traverse carbon atoms (the acyl chain should be mostly sp3/sp2 carbon backbone).
            if nbr.GetAtomicNum() == 6:
                if nbr_idx not in visited:
                    count += dfs_count_atoms(nbr_idx, visited)
        return count

    # In the acyl fragment, do not go through the carbonyl oxygen (which is part of the C=O).
    # Start DFS from the acyl carbon.
    visited = set()
    acyl_carbons = dfs_count_atoms(acyl_carbon_idx, visited)
    
    # For a short-chain fatty acid, we expect the acyl chain to have a limited number of carbons.
    # We set the threshold to 6 carbons (including the carbonyl carbon). Examples include
    # acetyl-, propionyl-, butyryl-, etc.
    if acyl_carbons > 6:
        return False, f"Acyl chain contains {acyl_carbons} carbons, too long for a short-chain fatty acid"
    if acyl_carbons < 2:
        return False, "Acyl chain too short (less than 2 carbons) to be a fatty acid"

    # Next, check for a CoA moiety.
    # A common fragment within CoA is the adenine ring. Many of our examples contain the pattern:
    # "n1cnc2c(N)ncnc12". (This is a simplified representation of the adenine moiety.)
    adenine_smarts = "n1cnc2c(n)cnc12"
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine fragment typical of CoA detected"
    
    # Additional optional checks could be done (e.g., verifying the presence of phosphates) 
    # but for this classification we assume the above criteria are sufficient.
    
    return True, f"Found thioester bond with an acyl chain of {acyl_carbons} carbons and a CoA moiety"

# (Optional) For testing, one might run examples:
if __name__ == "__main__":
    # Example: (R)-3-hydroxypentanoyl-CoA (one example provided)
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C[C@@H](CC)O)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
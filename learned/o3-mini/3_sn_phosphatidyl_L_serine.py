"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
  That is, a glycerol backbone with three –OH groups where the sn-3 hydroxyl is phosphorylated and then
  esterified with L-serine, and the sn-1 and sn-2 hydroxyls are esterified with fatty acyl chains.
  
This program uses RDKit to check that the given SMILES string has a phosphoserine head group 
and at least two acyl ester substituents.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    It verifies that the molecule contains:
      - A phosphoserine head group (a phosphate ester of L-serine attached via a glycerol fragment)
      - Two acyl ester substituents (as fatty acyl chains of at least 4 carbons) attached at the sn-1 and sn-2 positions.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the phosphoserine head group
    ps_pattern1 = Chem.MolFromSmarts("COP(=O)(O)OC[C@H](N)C(=O)O")
    ps_pattern2 = Chem.MolFromSmarts("COP(=O)(O)OC[C@@H](N)C(=O)O")
    
    # Check if either pattern is found
    ps_matches = mol.GetSubstructMatches(ps_pattern1)
    if not ps_matches:
        ps_matches = mol.GetSubstructMatches(ps_pattern2)
    if not ps_matches:
        return False, "Phosphoserine head group not found"
    
    # Record the atom indices for the phosphoserine head group from the first match (to avoid double counting)
    ps_atom_indices = set(ps_matches[0])
    
    # Define an ester SMARTS pattern that should capture the ester bond (carbonyl carbon and the ester oxygen)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Define a helper DFS function to calculate maximum contiguous carbon chain length starting from an atom.
    # We pass along a visited set to avoid loops.
    def get_chain_length(atom_idx, parent_idx, visited):
        """
        Recursively calculates the length of the longest contiguous chain (only considering carbons)
        starting from the given atom. 'visited' prevents cycles.
        """
        # Add the current atom to the visited set (for this branch)
        current_visited = visited.union({atom_idx})
        current_atom = mol.GetAtomWithIdx(atom_idx)
        max_length = 1  # Count the current atom
        for nb in current_atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == parent_idx or nb_idx in current_visited:
                continue
            if nb.GetAtomicNum() == 6: # Must be a carbon for an acyl chain
                branch_length = 1 + get_chain_length(nb_idx, atom_idx, current_visited)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length
    
    valid_acyl_count = 0
    # To avoid counting the same ester functionality twice, record carbonyl indices that have been processed.
    seen_carbonyls = set()
    
    # Process each ester match
    for match in ester_matches:
        # According to our SMARTS, match[0] is the carbonyl carbon and match[1] is the ester oxygen.
        carbonyl_idx = match[0]
        oxygen_idx = match[1]
        # Skip if this ester overlaps with the phosphoserine head group 
        if carbonyl_idx in ps_atom_indices or oxygen_idx in ps_atom_indices:
            continue
        
        # In the ester bond R-C(=O)-O-R, the acyl chain is tethered to the carbonyl carbon (excluding the oxygen).
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        acyl_start = None
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetIdx() == oxygen_idx:
                continue
            if nb.GetAtomicNum() == 6:  # likely the start of an acyl chain
                acyl_start = nb.GetIdx()
                break
        if acyl_start is None:
            continue  # no acyl chain found for this ester
        # Avoid double-counting the same carbonyl
        if carbonyl_idx in seen_carbonyls:
            continue
        seen_carbonyls.add(carbonyl_idx)
        
        # Calculate the maximum contiguous chain length using our DFS helper.
        chain_length = get_chain_length(acyl_start, carbonyl_idx, set())
        # We require at least 4 carbons to count as a valid fatty acyl chain.
        if chain_length >= 4:
            valid_acyl_count += 1

    if valid_acyl_count < 2:
        return False, f"Expected at least 2 acyl ester substituents on sn-1 and sn-2; found {valid_acyl_count}"
    
    # Optional: Check molecular weight as a rough filter (many phosphatidylserines typically have mol wt > 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a phosphatidylserine"
    
    return True, "Molecule contains a phosphoserine head group with two acyl ester substituents on a glycerol backbone"

# Optional test (uncomment for testing):
# test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCC"
# print(is_3_sn_phosphatidyl_L_serine(test_smiles))
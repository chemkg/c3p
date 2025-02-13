"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax Ester
Definition: A fatty acid ester resulting from the condensation of the carboxy group of a fatty acid 
with the alcoholic hydroxy group of a fatty alcohol.
"""

from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    Wax esters are characterized by a single ester linkage (-C(=O)-O-) connecting a fatty acid and a fatty alcohol.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a wax ester, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: perform a DFS (depth-first search) on the molecule starting from a given atom index,
    # avoiding traversal through any atom whose index is in the 'forbidden' set.
    def dfs_component(start_idx, forbidden):
        visited = set()
        stack = [start_idx]
        while stack:
            cur_idx = stack.pop()
            if cur_idx in visited:
                continue
            visited.add(cur_idx)
            atom = mol.GetAtomWithIdx(cur_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in visited and nbr_idx not in forbidden:
                    stack.append(nbr_idx)
        return visited

    candidate_ester = None  # Will store tuple: (ester_oxygen_idx, acid_carbon_idx, alcohol_carbon_idx)
    
    # Loop over all oxygen atoms to find candidate ester oxygens.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Look for oxygen that is connected to exactly two heavy atoms.
        if atom.GetDegree() != 2:
            continue
        
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
        
        # Identify candidate: one neighbor (acid candidate) must be a carbon that is double-bonded to an oxygen.
        acid_candidate = None
        alcohol_candidate = None
        for nbr in neighbors:
            if nbr.GetAtomicNum() != 6:
                continue
            # Check if this carbon is attached via a double bond to an oxygen (other than our candidate oxygen).
            found_doubleO = False
            for bond in nbr.GetBonds():
                # Get the other atom in this bond.
                other = bond.GetOtherAtom(nbr)
                if other.GetIdx() == atom.GetIdx():
                    continue
                if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    found_doubleO = True
                    break
            if found_doubleO:
                acid_candidate = nbr.GetIdx()
            else:
                alcohol_candidate = nbr.GetIdx()
        
        if acid_candidate is not None and alcohol_candidate is not None:
            if candidate_ester is not None:
                # More than one candidate ester group found.
                return False, "More than one ester group found; not a simple wax ester"
            candidate_ester = (atom.GetIdx(), acid_candidate, alcohol_candidate)
    if candidate_ester is None:
        return False, "No ester group found that meets the wax ester pattern"
    
    ester_oxygen_idx, acid_carbon_idx, alcohol_carbon_idx = candidate_ester

    # Define a minimum threshold for the number of carbon atoms in each chain.
    MIN_CARBONS = 6
    
    # For the fatty alcohol part: traverse starting from alcohol_carbon_idx, 
    # blocking through the bridging ester oxygen.
    alcohol_fragment = dfs_component(alcohol_carbon_idx, forbidden={ester_oxygen_idx})
    alcohol_carbons = [idx for idx in alcohol_fragment if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    
    if len(alcohol_carbons) < MIN_CARBONS:
        return False, f"Fatty alcohol portion too short; found {len(alcohol_carbons)} carbon(s)"
    
    # For the fatty acid part: traverse starting from the carbonyl carbon (acid_carbon_idx),
    # blocking the bridging ester oxygen.
    acid_fragment = dfs_component(acid_carbon_idx, forbidden={ester_oxygen_idx})
    acid_carbons = [idx for idx in acid_fragment if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    
    if len(acid_carbons) < MIN_CARBONS:
        return False, f"Fatty acid portion too short; found {len(acid_carbons)} carbon(s)"
    
    # If our molecule contains additional ester groups (other than the one we found),
    # it may not be a simple wax ester. We can quickly check for any other ester pattern.
    # We use a simple SMARTS for an ester functional group (excluding the one we already identified).
    ester_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ester_matches = mol.GetSubstructMatches(ester_smarts)
    if len(ester_matches) > 1:
        return False, "More than one ester group detected in the molecule"

    return True, (f"Wax ester detected: fatty acid chain with {len(acid_carbons)} C(s) and fatty alcohol chain with "
                  f"{len(alcohol_carbons)} C(s) connected via a single ester linkage.")

# Example usage:
if __name__ == "__main__":
    test_smiles = "O(CCCCCCCC/C=C\\CCCCCC)C(=O)CCCCCCCCCCCCCCCCCCC"  # Palmitoleyl arachidate
    result, message = is_wax_ester(test_smiles)
    print(result, message)
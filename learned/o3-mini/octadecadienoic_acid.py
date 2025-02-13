"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
#!/usr/bin/env python
"""
Classifies: Octadecadienoic acid
Definition: Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.
In this implementation we require that:
  - The molecule contains a carboxylic acid group (C(=O)[OH])
  - Starting from the carboxyl carbon (the C in the acid group) there is a single, linear, unbranched chain of 18 carbon atoms.
  - Counting only the bonds linking consecutive carbons in that chain (and ignoring the C=O of the acid group),
    exactly 2 of these bonds are C=C double bonds.
Note: This procedure is a simplified approach and may not cover all edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    Octadecadienoic acid is defined as any straight-chain, C18 fatty acid with exactly two C=C double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as octadecadienoic acid, False otherwise.
        str: A textual reason for the decision.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find the carboxylic acid group using SMARTS: C(=O)[OH]
    acid_smarts = "C(=O)[OH]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
        
    # Assume the first match is our acid group.
    # Identify the acid carbon (the one with atomic number 6) from the match.
    acid_carbon = None
    for atom_idx in acid_matches[0]:
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6:
            acid_carbon = atom_idx
            break
    if acid_carbon is None:
        return False, "Carboxylic acid pattern found but no acid carbon identified"
    
    # Now we try to locate the straight chain of carbons starting from the acid carbon.
    # In a straight chain, the acid carbon should be attached to exactly one carbon (the chain continuation).
    chain = [acid_carbon]  # this list of atom indices will represent our candidate fatty acid chain
    acid_atom = mol.GetAtomWithIdx(acid_carbon)
    chain_neighbors = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(chain_neighbors) != 1:
        return False, f"Acid carbon has {len(chain_neighbors)} carbon neighbors, expected exactly one for a straight-chain fatty acid"
    next_atom = chain_neighbors[0]
    chain.append(next_atom)
    
    # Walk down the chain in a linear, unbranched fashion.
    # At every step, from the current atom, choose the next carbon neighbor not already in our chain.
    while len(chain) < 18:
        current_atom = mol.GetAtomWithIdx(chain[-1])
        # Get neighbors that are carbon and not already in the chain.
        cand_neighbors = [nbr.GetIdx() for nbr in current_atom.GetNeighbors() 
                          if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in chain]
        if len(cand_neighbors) == 0:
            break  # chain ended prematurely
        if len(cand_neighbors) > 1:
            return False, f"Branching detected at atom index {chain[-1]} (chain is not straight)"
        chain.append(cand_neighbors[0])
    
    if len(chain) != 18:
        return False, f"Chain length is {len(chain)}, expected 18 carbons for an octadecadienoic acid"
    
    # Verify that the terminal carbon is a methyl group.
    terminal_atom = mol.GetAtomWithIdx(chain[-1])
    # For a methyl group (in a straight chain), the only heavy-atom neighbor should be the preceding chain carbon.
    terminal_carbons = [nbr.GetIdx() for nbr in terminal_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(terminal_carbons) != 1:
        return False, "Terminal carbon is not a methyl group (unexpected extra carbon neighbors)"
    
    # Count the carbon–carbon double bonds along the chain.
    # We ignore the carbonyl double bond in the acid group.
    dbl_bond_count = 0
    for i in range(len(chain) - 1):
        bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
        if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
            dbl_bond_count += 1
    if dbl_bond_count != 2:
        return False, f"Found {dbl_bond_count} carbon–carbon double bond(s) in the chain, expected exactly 2"
    
    return True, "Molecule is a straight-chain C18 fatty acid with 2 C=C bonds"

# Example usage:
if __name__ == "__main__":
    # Try one of the given examples:
    test_smiles = "CCCCCC\\C=C\\C=C/CCCCCCCC(O)=O"  # 9-cis,11-trans-octadecadienoic acid
    result, reason = is_octadecadienoic_acid(test_smiles)
    print(result, reason)
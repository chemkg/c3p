"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid containing at least one C=C double bond.
A valid olefinic fatty acid (for this classifier) is a long‐chain (>=12 C atoms) acyclic molecule that:
  - Contains exactly one terminal carboxylic acid group (C(=O)[O;H,-] with the acid C having exactly one carbon neighbor)
  - Contains at least one C=C double bond.
  - The majority (>70%) of its non‐H atoms are carbons, reflecting that it is essentially an unadorned fatty acyl chain.
  
Note: This is a heuristic and may fail for edge cases.
"""

from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    Olefinic fatty acids are defined (here) as acyclic long‐chain carboxylic acids 
    that have a terminal carboxylic acid group (the acid C is connected to only one carbon), 
    contain at least one C=C double bond, and where the majority of non‐hydrogen atoms are carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if the molecule contains any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a simple acyclic fatty acid"
    
    # Check for a carboxylic acid group.
    # The SMARTS looks for a carbon double-bonded to O and single-bonded to an -OH (or deprotonated O).
    acid_smarts = "C(=O)[O;H1,-]"  # the first atom is the acid carbon.
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # We require exactly one terminal carboxylic acid group.
    terminal_acid_found = False
    for match in acid_matches:
        acid_c_idx = match[0]  # first atom: the carbonyl carbon of the acid.
        acid_atom = mol.GetAtomWithIdx(acid_c_idx)
        # Count carbon neighbors (ignoring oxygens); a terminal acid carbon should have exactly one carbon neighbor.
        carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            break
    if not terminal_acid_found:
        return False, "Carboxylic acid group is not terminal; not a typical fatty acid"
    
    # Count total number of carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 12:
        return False, f"Too few carbon atoms ({len(c_atoms)} found, need at least 12) to be a long-chain fatty acid"
    
    # Check for the presence of at least one carbon-carbon double bond.
    # This SMARTS matches a C=C bond (and will not confuse it with a C=O).
    olefin_smarts = "C=C"
    olefin_pattern = Chem.MolFromSmarts(olefin_smarts)
    if not mol.HasSubstructMatch(olefin_pattern):
        return False, "No C=C double bonds found; not an olefinic fatty acid"
    
    # Heuristic: Most non-hydrogen atoms in a fatty acid should be carbons.
    non_H_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_nonH = len(non_H_atoms)
    n_carbons = sum(1 for atom in non_H_atoms if atom.GetAtomicNum() == 6)
    if n_nonH == 0 or (n_carbons / n_nonH) < 0.7:
        return False, "The molecule does not appear to be a simple fatty acyl chain (low fraction of carbon atoms)"
    
    return True, "Contains a terminal carboxylic acid group, is acyclic and long-chained (>=12 C atoms), and has at least one C=C double bond, classifying it as an olefinic fatty acid"

# Example usage:
if __name__ == '__main__':
    test_smiles = "CCCCCC\\C=C/CC\\C=C/CCCC(O)=O"  # (example: a simplified olefinic fatty acid)
    result, reason = is_olefinic_fatty_acid(test_smiles)
    print(result, reason)
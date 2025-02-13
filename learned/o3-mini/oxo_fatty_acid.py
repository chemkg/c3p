"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo Fatty Acid
Definition: Any fatty acid containing at least one additional aldehydic or ketonic (oxo) group 
in addition to a terminal free carboxylic acid group. Additionally, the molecule should exhibit 
a predominantly acyclic (or contiguous) alkyl chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid in this context is defined as a (possibly short‐ or long‐chain) fatty acid 
    with a terminal (free) carboxylic acid group and at least one extra oxo group (ketone or aldehyde)
    not belonging to the carboxylic acid functionality. In addition the contiguous acyclic chain emanating 
    from the acid should be a large fraction of the molecule’s carbons.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if molecule meets oxo fatty acid criteria, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (A) Look for a terminal free carboxylic acid.
    # Standard SMARTS for carboxylic acid: [CX3](=O)[OX2H1]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group detected; not a fatty acid"
    
    # Identify terminal acid groups.
    # We require that the acid carbon (first atom in SMARTS) is terminal, meaning it is bonded to exactly one carbon.
    terminal_acid_indices = []  # store acid carbon indices that qualify
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_atom = mol.GetAtomWithIdx(acid_carbon_idx)
        # Count neighboring carbon atoms (ignore oxygens).
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_indices.append((acid_carbon_idx, carbon_neighbors[0]))
    
    if not terminal_acid_indices:
        return False, "Carboxylic acid group is not terminal; fatty acids require a free (terminal) acid group"
    
    # For further analysis, choose the first terminal acid. Its unique neighbor (acyl chain start) is:
    acid_carbon_idx, chain_start_idx = terminal_acid_indices[0]
    
    # (B) Check that there is at least one additional oxo group (aldehyde or ketone) outside the acid.
    # SMARTS for ketone: a carbonyl with carbons on both sides.
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    # SMARTS for aldehyde: a carbonyl with at least one hydrogen attached to the carbonyl carbon.
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H](=O)")
    
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    additional_oxo_found = False
    # For any ketone, ensure that the carbonyl carbon is not our terminal acid carbon.
    for match in ketone_matches:
        # In our ketone pattern the 2nd atom is the carbonyl carbon.
        if match[1] != acid_carbon_idx:
            additional_oxo_found = True
            break
    if not additional_oxo_found:
        for match in aldehyde_matches:
            if match[1] != acid_carbon_idx:
                additional_oxo_found = True
                break
    if not additional_oxo_found:
        return False, "No additional oxo (aldehyde/ketone) group detected outside the acid"

    # (C) Examine the acyclic alkyl chain.
    # Count total carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 6:
        return False, "Too few carbon atoms to be considered a fatty acid"
    
    # Now, from the chain starting point (neighboring the acid carbon), extract the contiguous acyclic carbon chain.
    # We restrict our DFS to carbon atoms that are not in any ring.
    acyclic_chain_set = set()
    stack = [chain_start_idx]
    
    while stack:
        current_idx = stack.pop()
        if current_idx in acyclic_chain_set:
            continue
        atom = mol.GetAtomWithIdx(current_idx)
        # Only consider carbon atoms not in a ring.
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            acyclic_chain_set.add(current_idx)
            # Add neighbors that are carbons (if not the acid carbon which we already used).
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() not in acyclic_chain_set):
                    # Continue only if neighbor is non-ring.
                    if not nbr.IsInRing():
                        stack.append(nbr)
    chain_length = len(acyclic_chain_set)
    
    # Many fatty acids have a long contiguous chain. 
    # To allow even medium-chain oxo fatty acids we use a moderate lower bound (~4 carbons in the chain)
    # and require that the chain comprises at least 50% of all carbons.
    if chain_length < 4:
        return False, f"Contiguous acyclic chain is too short ({chain_length} carbons found)"
    if chain_length / total_carbons < 0.5:
        return False, ("The acyclic chain only represents {:.0%} of all carbons; a fatty acid should be predominantly "
                       "aliphatic".format(chain_length / total_carbons))
    
    return True, ("Contains a terminal free carboxylic acid group, an additional oxo (aldehyde/ketone) group "
                  "outside that acid, and a sufficiently long acyclic carbon chain ({} of {} carbons)".format(chain_length, total_carbons))


# Example usage:
if __name__ == "__main__":
    # Test example: (6E,8Z)-5-oxooctadecadienoic acid
    test_smiles = "O=C(CCCC(O)=O)/C=C/C=C\\CCCCCCCCC"
    result, reason = is_oxo_fatty_acid(test_smiles)
    print("Test molecule classification:", result)
    print("Reason:", reason)
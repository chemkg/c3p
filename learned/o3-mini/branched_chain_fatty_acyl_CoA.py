"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any branched-chain fatty acid.
The improved algorithm:
  1. Checks for presence of a CoA moiety (via an adenine-derived substructure).
  2. Checks for a thioester bond linking the fatty acyl chain to CoA.
  3. Starting from the carbonyl’s fatty acyl neighbor (i.e. not the S partner),
     it traverses only carbon atoms (excluding the carbonyl) to collect the acyl chain.
  4. Requires that the chain is long enough and that at least one carbon shows branching 
     (i.e. it has more than two carbon neighbors within the chain).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    The improved algorithm:
      1. Parses the SMILES.
      2. Requires an adenine substructure from CoA.
      3. Requires a thioester bond [CX3](=O)[SX2] linking the fatty acyl chain.
      4. From each thioester, identifies the fatty acyl branch (obtained from the carbonyl atom's neighbor
         that is not the sulfur) and performs a DFS (restricted to carbon atoms and not going back to the carbonyl)
         to capture the fatty acyl chain.
      5. Checks that the chain is long (at least 3 atoms) and that at least one carbon in the chain has more than 
         2 neighbors (within the chain) indicating branching.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a branched-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification.
    """
    
    # Step 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2. Check for the presence of a CoA moiety.
    # Searching for the adenine substructure found in Coenzyme A.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"
    
    # Step 3. Check for the thioester linkage (carbonyl carbon linked to sulfur).
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond not found"

    # Process each thioester match until one yields a branched fatty acyl chain.
    for match in thioester_matches:
        # match: first atom is the carbonyl C, second is the S connected to CoA.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Step 4. Identify the fatty acyl neighbor: a neighbor of the carbonyl that is not the sulfur
        fatty_acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue  # Skip the S that goes to CoA.
            # Skip the double-bonded oxygen of the carbonyl.
            if nbr.GetAtomicNum() == 8:
                continue
            # We assume the fatty acyl chain begins at this neighbor.
            fatty_acyl_neighbor = nbr
            break
        
        if fatty_acyl_neighbor is None:
            continue  # Try next thioester match if fatty acyl side is not found.
        
        # Step 5. Build the fatty acyl chain by traversing only carbon atoms.
        # We start from the fatty acyl neighbor and do not step back to the carbonyl.
        acyl_chain = set()
        stack = [fatty_acyl_neighbor.GetIdx()]
        while stack:
            current_idx = stack.pop()
            if current_idx in acyl_chain:
                continue
            atom = mol.GetAtomWithIdx(current_idx)
            # Only consider carbon atoms (atomic number 6).
            if atom.GetAtomicNum() != 6:
                continue
            acyl_chain.add(current_idx)
            for nbr in atom.GetNeighbors():
                # Do not go back to the carbonyl atom.
                if nbr.GetIdx() == carbonyl_idx:
                    continue
                # Continue only if neighbor is carbon.
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in acyl_chain:
                    stack.append(nbr.GetIdx())
        
        # Require a minimum chain length (e.g., at least 3 carbons) to qualify as a fatty acyl chain.
        if len(acyl_chain) < 3:
            continue  # Try next match
        
        # Step 6. Check for branching within the acyl chain.
        # In a linear chain, internal carbons have exactly 2 carbon neighbors (within the chain) and terminal have 1.
        branch_found = False
        for idx in acyl_chain:
            atom = mol.GetAtomWithIdx(idx)
            # Count only neighbors that are in the acyl_chain.
            n_chain_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in acyl_chain)
            if n_chain_neighbors > 2:
                branch_found = True
                break
        
        if branch_found:
            return True, "Branched-chain fatty acyl-CoA detected"
        else:
            # If no branch is found in this acyl chain, it is most likely linear.
            return False, "Acyl chain appears linear; no branch detected"
    
    # If none of the thioester groups yield a valid fatty acyl chain.
    return False, "No fatty acyl chain meeting criteria found"


# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES strings.
    sample_smiles = "CC(C)CCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(sample_smiles)
    print(result, reason)
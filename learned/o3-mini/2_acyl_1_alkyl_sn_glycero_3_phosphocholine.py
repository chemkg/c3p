"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine

Definition:
An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl groups are located 
at positions 1 and 2 respectively. These molecules consist of three parts:
    • A phosphocholine headgroup (characterized by the pattern OP(=O)([O-])OCC[N+](C)(C)C)
    • One acyl ester chain (O–C(=O)… motif, exactly one must be present)
    • One O-alkyl (ether) chain that is not part of the ester substructure and directly connected to an oxygen
       that leads to a long uninterrupted aliphatic chain.

This program uses RDKit substructure searches and a recursive chain length search to distinguish 
this lipid type from other similar lipids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_alkyl_chain(atom, visited=None):
    """
    Recursively finds the longest chain (in terms of the number of carbon atoms) 
    starting from the given carbon atom. Only aliphatic, sp3 carbons (atomic number 6) 
    that are not part of a carbonyl (i.e. not double-bonded to an oxygen) are traversed.
    
    Args:
        atom: RDKit Atom to start the search.
        visited: Set of atom indices already visited.
    
    Returns:
        int: Length (number of bonds) of the longest alkyl chain found.
    """
    if visited is None:
        visited = set()
    visited.add(atom.GetIdx())
    max_length = 0
    for nbr in atom.GetNeighbors():
        # Continue only for carbon atoms not yet visited.
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
            # Check if the carbon is part of a carbonyl (skip if so).
            is_carbonyl = False
            for bond in nbr.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8:
                        is_carbonyl = True
                        break
            if is_carbonyl:
                continue
            branch_length = 1 + longest_alkyl_chain(nbr, visited.copy())
            if branch_length > max_length:
                max_length = branch_length
    return max_length

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Improved heuristic logic:
      1. Confirm the presence of the phosphocholine headgroup. We search for the pattern 
         OP(=O)([O-])OCC[N+](C)(C)C.
      2. Count the number of ester groups by finding substructures matching O[C](=O); exactly one must be present 
         to define the acyl chain.
      3. Identify an O-alkyl (ether) linkage that is not part of the ester substructure and is not connected to phosphorus.
         Then, for each candidate oxygen, check that one of its carbon neighbors leads to a long alkyl chain (chain length >= 8).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine.
        str: A reason string for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the phosphocholine headgroup.
    phos_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phos_pattern):
        return False, "Missing phosphocholine headgroup (pattern OP(=O)([O-])OCC[N+](C)(C)C not found)"
    
    # 2. Look for the acyl ester group (O–C(=O) motif).
    ester_pattern = Chem.MolFromSmarts("O[C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "Missing acyl ester group (pattern O[C](=O) not found)"
    if len(ester_matches) != 1:
        return False, f"Expected exactly one acyl ester group but found {len(ester_matches)}; may be a diacyl lipid."
    
    # Record oxygen atoms that belong to the ester group. 
    # In our pattern, the first atom (index 0) is the oxygen attached to the acyl chain.
    ester_oxygen_indices = set(match[0] for match in ester_matches)
    
    # 3. Identify an O-alkyl (ether) linkage.
    ether_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:  # Only consider oxygen atoms.
            continue
        # Exclude oxygens that are part of the ester substructure.
        if atom.GetIdx() in ester_oxygen_indices:
            continue
        # Exclude oxygens connected directly to phosphorus (as in the phosphate headgroup).
        if any(neigh.GetAtomicNum() == 15 for neigh in atom.GetNeighbors()):
            continue
        
        # For each remaining oxygen, examine neighboring carbon atoms.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            # Use the recursive function to compute the longest chain from this carbon.
            chain_length = longest_alkyl_chain(nbr)
            if chain_length >= 8:
                ether_found = True
                break
        if ether_found:
            break
    
    if not ether_found:
        return False, "Missing or too short O-alkyl (ether) chain; no oxygen found leading to an alkyl chain of sufficient length."
    
    return True, ("Contains a phosphocholine headgroup, exactly one acyl ester group, "
                   "and one O-alkyl (ether) chain with a long uninterrupted alkyl chain. "
                   "This is consistent with 2-acyl-1-alkyl-sn-glycero-3-phosphocholine.")

# For testing purposes:
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCC"
    result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(test_smiles)
    print(f"Test result: {result}\nReason: {reason}")
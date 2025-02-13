"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate ester 
Definition: Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).

Improvement rationale:
  • Instead of merely checking for the substring "CCCCCCCC(=O)O", we search for ester groups.
  • For each ester group we identify the acyl fragment attached to the carbonyl carbon.
  • We then traverse that chain (ensuring it is an unbranched aliphatic chain with only single bonds)
    and count the number of carbons. An octanoate ester must have an acyl group with exactly eight carbons
    (the carbonyl carbon counts as one).
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def _count_linear_acyl_chain(mol: rdchem.Mol, carbonyl, acyl_neighbor):
    """
    Starting from the carbonyl carbon (already counted as 1) and its acyl side neighbor,
    follow the chain along single bonds between sp3 carbons that form a linear (nonbranching) alkyl chain.
    Return the total number of carbons in the acyl chain (including the carbonyl carbon).
    If a branch is encountered, then the chain is not linear and we return None.
    """
    count = 1  # starting with the carbonyl carbon
    prev_atom = carbonyl
    current = acyl_neighbor
    # Check that the acyl neighbor is a saturated carbon
    if current.GetAtomicNum() != 6:
        return None
    # Increase count for the acyl neighbor
    count += 1

    while True:
        # Find carbon neighbors of current that are not the one we came from.
        next_carbons = []
        for nbr in current.GetNeighbors():
            if nbr.GetIdx() == prev_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                # Also ensure the bond is a single bond
                bond = mol.GetBondBetweenAtoms(current.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == rdchem.BondType.SINGLE:
                    next_carbons.append(nbr)
        if len(next_carbons) == 0:
            # End of chain reached
            break
        elif len(next_carbons) == 1:
            # Linear continuation
            prev_atom = current
            current = next_carbons[0]
            count += 1
        else:
            # Branch detected; not a linear fatty acid chain.
            return None
    return count

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is defined as any fatty acid ester in which the acid component is octanoic acid.
    Octanoic acid contains 8 carbons (CH3(CH2)6COOH); hence the acyl group has exactly eight carbons,
    counting the carbonyl carbon.
    
    The algorithm works as follows:
      1. Parse the SMILES to obtain the molecule.
      2. Find ester groups using the SMARTS pattern for an ester fragment [#6](=O)[O].
      3. For each match (the carbonyl carbon and the ester oxygen), identify the acyl side by taking the
         carbonyl neighbor that is not the ester oxygen.
      4. Traverse this acyl chain as a linear series of sp3 carbons (single bonds only) and count the total
         number of carbons. If the count equals 8, then the acid component must be octanoic acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an octanoate ester moiety, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for an ester group:
    # The pattern [#6](=O)[O] matches a carbonyl carbon (atomic number 6, with three connections)
    # attached via a double bond to oxygen and single bond to an oxygen.
    ester_smarts = "[#6X3](=O)[OX2H0]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return None, None  # something went wrong
    
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No ester group found in molecule"
    
    # Look at each ester match and test the acyl chain length.
    for match in matches:
        # In the SMARTS, match[0] is the carbonyl carbon and match[1] is the ester oxygen.
        carbonyl = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[1])
        
        # Identify the acyl-chain neighbor of the carbonyl carbon (ignore the oxygen).
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() if nbr.GetIdx() != ester_oxygen.GetIdx() and nbr.GetAtomicNum() == 6]
        if not acyl_neighbors:
            continue  # no acyl side found in this ester
        
        # We assume there is only one acyl side (if more, each will be checked)
        acyl_neighbor = acyl_neighbors[0]
        chain_len = _count_linear_acyl_chain(mol, carbonyl, acyl_neighbor)
        if chain_len is None:
            # chain is not linear so skip
            continue
        
        # For octanoate ester the acyl chain length must be exactly 8 carbons (including the carbonyl carbon).
        if chain_len == 8:
            return True, "Molecule contains an octanoate ester moiety"
    
    return False, "Molecule does not contain an octanoate ester moiety"

# Example usage:
if __name__ == "__main__":
    test_smiles_list = [
        "CCCCCCCC(=O)OCC",  # ethyl octanoate
        "CCCCCCCC(=O)OC",   # methyl octanoate
        "CC(=O)OCCCCCCCC(=O)O",  # a molecule with one octanoate ester and one acetyl group
        "CCCCCCCC(=O)OC[C@H](O)CO", # 1-octanoyl-sn-glycerol
        # Some false positive example (should not match):
        "CCC(=O)O[C@H]1[C@H](C)O[C@H](C[C@@]1(C)O)...",  # truncated leucomycin like SMILES
    ]
    for smi in test_smiles_list:
        result, reason = is_octanoate_ester(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")
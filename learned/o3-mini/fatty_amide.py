"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: Fatty Amide
Definition: A monocarboxylic acid amide derived from a fatty acid.
A fatty amide is recognized by an amide bond (C(=O)N) where the acyl (carbonyl) portion is derived from a fatty acid
(i.e. contains a long aliphatic chain of at least 8 carbon atoms) while the amine portion is relatively small (no long aliphatic chain).
"""

from rdkit import Chem

def longest_aliphatic_chain_length(atom, visited):
    """
    Recursively computes the length (number of carbon atoms) of the longest _linear_ (acyclic, non-aromatic) chain
    starting from the given atom. We only traverse carbon atoms that are not aromatic and not in a ring.
    'visited' is a set of atom indices to prevent cycles.
    """
    max_length = 1  # count current atom
    for nbr in atom.GetNeighbors():
        if (nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() and not nbr.IsInRing() 
            and nbr.GetIdx() not in visited):
            # Copy visited so that different branches do not interfere.
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            branch_length = 1 + longest_aliphatic_chain_length(nbr, new_visited)
            if branch_length > max_length:
                max_length = branch_length
    return max_length

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide.

    A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
    That means that:
      (1) an amide bond (C(=O)N) is present,
      (2) the acyl (C=O side) must come from a fatty acid and thus has an aliphatic chain of at least 8 carbon atoms (including the carbonyl carbon),
      (3) at least one of the substituents on the amide nitrogen is not fatty (i.e. does not have a long aliphatic chain – here we require all chains off the nitrogen be shorter than 8 carbons).
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule meets the fatty amide criteria, False otherwise.
        str: Message giving the reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for a simple amide group: the pattern returns carbonyl carbon, oxygen and nitrogen.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if amide_smarts is None:
        return False, "Error generating amide SMARTS pattern"
        
    matches = mol.GetSubstructMatches(amide_smarts)
    if not matches:
        return False, "No amide group found in the molecule"

    # Define thresholds: fatty acyl chain must be at least 8 carbons.
    fatty_chain_threshold = 8
    # And the substituent on the nitrogen must be short; here we require less than 8 carbons.
    amine_chain_max = 7

    # Process each amide match
    for match in matches:
        # match: (carbonyl carbon, oxygen, amide nitrogen)
        carbonyl = mol.GetAtomWithIdx(match[0])
        oxy = mol.GetAtomWithIdx(match[1])
        amide_nitrogen = mol.GetAtomWithIdx(match[2])
        
        # From the carbonyl carbon, get neighbors that represent the fatty acyl chain.
        # Exclude the carbonyl oxygen and amide nitrogen.
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors()
                          if nbr.GetIdx() not in {oxy.GetIdx(), amide_nitrogen.GetIdx()}
                          and nbr.GetAtomicNum() == 6
                          and not nbr.GetIsAromatic()
                          and not nbr.IsInRing()]
        if not acyl_neighbors:
            # Cannot find a carbon chain on the acyl side.
            continue

        # Compute maximum chain length from carbonyl: include the carbonyl carbon itself.
        acyl_chain_length = 0
        for nbr in acyl_neighbors:
            chain_length = 1 + longest_aliphatic_chain_length(nbr, visited={carbonyl.GetIdx()})
            if chain_length > acyl_chain_length:
                acyl_chain_length = chain_length

        if acyl_chain_length < fatty_chain_threshold:
            # This amide does not qualify as having a sufficiently long fatty acyl chain.
            continue

        # Now check the amide nitrogen side: ensure that no substituent on the nitrogen (other than the carbonyl)
        # is too long (i.e. fatty). We want the amine part to be relatively small.
        valid_amine = True
        for nbr in amide_nitrogen.GetNeighbors():
            if nbr.GetIdx() == carbonyl.GetIdx():
                continue  # skip the carbonyl connection
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() and not nbr.IsInRing():
                chain_len = 1 + longest_aliphatic_chain_length(nbr, visited={amide_nitrogen.GetIdx()})
                if chain_len >= fatty_chain_threshold:
                    valid_amine = False
                    break
        if not valid_amine:
            # This amide's nitrogen bears a long aliphatic chain that makes it likely part of a lipid (e.g. ceramide)
            # and not a typical fatty amide where only the acyl group is long.
            continue

        # If we reach here, we found an amide bond whose acyl chain is long enough and whose nitrogen substituents are not fatty.
        return True, f"Found fatty amide: amide group with acyl chain length {acyl_chain_length} carbons and modest amine substituents."
        
    return False, "Amide group found but no qualifying fatty acyl chain or amine substituent criteria not met."
    
# Example usage:
if __name__ == "__main__":
    examples = [
        "CC(C)CC(=O)NCCc1c[nH]cn1",  # Dolichotheline (should be false: acyl chain too short)
        "CCCCCCCCC(=O)NCc1ccc(O)c(OC)c1",  # nonivamide (acyl chain length 9 --> true)
        "O=C(NCC(O)=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC"  # Docosahexaenoyl glycine (22 carbons acyl)
    ]
    for smi in examples:
        result, reason = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
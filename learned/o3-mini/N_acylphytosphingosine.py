"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
The code uses a SMARTS pattern able to capture the motif:
    acyl carbonyl – C(=O) attached to N
    which is attached to a sphingosine head of the form N–C(CO)–C(O)–C(O)
In addition the acyl fragment and the sphingoid chain (beyond the head) are required to be sufficiently long.
Note: The classification is based on heuristic substructure patterns.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    A N-acylphytosphingosine is defined as a ceramide where phytosphingosine is acylated at the nitrogen.
    
    The expected substructure is roughly:
       –[R]-C(=O)–N–C(CO)–C(O)–C(O)–[tail]
    where the acyl (fatty acyl) group (to the left of C(=O)) and the sphingoid tail
    (to the right of the three–carbon head) are long aliphatic chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is N-acylphytosphingosine, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for the core. This pattern looks for:
    #   [C](=O)[N][C](CO)[C](O)[C](O)
    # It encodes an acyl carbonyl (which should have a double bond O) bound to an amide N,
    # then the phosphingosine head made of three carbons, where the first carries "CO" (a CH2OH),
    # followed by two carbons with OH groups.
    core_smarts = "[C](=O)[N][C](CO)[C](O)[C](O)"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        # In case SMARTS did not compile (should not happen)
        return None, None

    matches = mol.GetSubstructMatches(core_pattern)
    if not matches:
        return False, "Core motif [C](=O)[N][C](CO)[C](O)[C](O) not found: not an N-acylphytosphingosine"

    # At this point, we have one or more candidate matches.
    # We require that at least one match passes further tests.
    for match in matches:
        # match is a tuple of 5 atom indices corresponding to:
        #   idx0: acyl carbon (carbonyl carbon)
        #   idx1: amide nitrogen
        #   idx2: first carbon of sphingosine head (should have a CH2OH substituent)
        #   idx3: second carbon (with OH)
        #   idx4: third carbon (with OH) that typically connects to the sphingoid tail.
        acylC = mol.GetAtomWithIdx(match[0])
        amideN = mol.GetAtomWithIdx(match[1])
        head1 = mol.GetAtomWithIdx(match[2])
        head2 = mol.GetAtomWithIdx(match[3])
        head3 = mol.GetAtomWithIdx(match[4])
        
        # Verify that the acyl carbon (idx0) has an oxygen connected by a double bond.
        # That oxygen should be the carbonyl oxygen.
        has_carbonyl_ox = False
        for nb in acylC.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(acylC.GetIdx(), nb.GetIdx())
            if nb.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_carbonyl_ox = True
                break
        if not has_carbonyl_ox:
            continue  # try next match
        
        # Check that the acyl group is sufficiently "fatty" (e.g., at least 6 carbons in the acyl chain)
        # We do a DFS from the acyl carbon. For this purpose we do not go into the amide (amideN) because
        # that leads into the sphingosine core.
        def dfs_count_carbons(atom, visited):
            count = 0
            if atom.GetAtomicNum() == 6:
                count = 1
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in visited:
                    continue
                # Do not go into the amide nitrogen (or its nearby atoms that are part of the core)
                if nb.GetIdx() in core_match_set:
                    continue
                # Only traverse through bonds that are single; this restricts to the aliphatic chain.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                visited.add(nb.GetIdx())
                count += dfs_count_carbons(nb, visited)
            return count

        # Create a set of indices that belong to the core match so we do not traverse back:
        core_match_set = set(match)
        # For acyl chain, start from the acyl carbon and follow neighbors except the amide nitrogen.
        acyl_chain_count = 0
        for nb in acylC.GetNeighbors():
            if nb.GetIdx() == amideN.GetIdx():
                continue
            # Traverse from this neighbor
            visited = set([acylC.GetIdx()])
            acyl_chain_count = max(acyl_chain_count, dfs_count_carbons(nb, visited))
        if acyl_chain_count < 6:
            # too short to be a fatty acyl chain
            continue  # try next match

        # Now check that the sphingoid tail (the chain attached to the third head carbon) is sufficiently long.
        # Note that head3 (match[4]) is expected to have an alkyl tail.
        tail_chain_count = 0
        for nb in head3.GetNeighbors():
            # Avoid going back to the core head (head2) by checking membership in the core
            if nb.GetIdx() in core_match_set:
                continue
            visited = set([head3.GetIdx()])
            tail_chain_count = max(tail_chain_count, dfs_count_carbons(nb, visited))
        if tail_chain_count < 6:
            continue  # too short sphingoid tail

        # If we have reached this point for one candidate match then we believe the structure to be N-acylphytosphingosine.
        return True, "Molecule contains an acylated phytosphingosine core with sufficiently long fatty acyl and sphingoid chains"
    
    # If none of the candidate matches passed the further checks, then return false
    return False, "Matches found for the core motif but did not meet acyl or sphingoid chain criteria"

# For informal testing uncomment below lines:
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
# print(is_N_acylphytosphingosine(test_smiles))
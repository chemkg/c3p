"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
The code first finds a core motif roughly corresponding to:
    –R–C(=O)–N–C(CO)[C](O)[C](O)–
It then verifies that:
  1. The acyl carbon (attached to the carbonyl) indeed has a double-bonded oxygen.
  2. Each of the three head carbons (the ones coming off the N) bear hydroxyl groups.
     The hydroxyl check is enhanced to detect both direct –OH attachment or if a CH2 group
     leads to an –OH.
  3. Both the fatty acyl chain (R–) and the sphingoid tail (attached via the third head carbon)
     are sufficiently long (at least 8 carbon atoms).
If any of these criteria fails then the molecule is rejected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    
    A N-acylphytosphingosine is (roughly) defined as a ceramide with a phytosphingosine base,
    having the following core pattern:
         R–C(=O)–N–C(CO)[C](O)[C](O)
    In addition, both the fatty acyl chain (R–) and the sphingoid tail (attached at the third head carbon)
    are required to be long (here at least 8 carbon atoms long) and each head carbon must carry a hydroxyl group.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if molecule is classified as an N-acylphytosphingosine, False otherwise.
       str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core SMARTS (ignoring chirality).
    # Note: [C](CO) will often match [C@@H](CO) as chirality is not enforced here.
    core_smarts = "[C](=O)[N][C](CO)[C](O)[C](O)"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return None, None  # Should not happen
    
    matches = mol.GetSubstructMatches(core_pattern)
    if not matches:
        return False, f"Core motif {core_smarts} not found: not an N-acylphytosphingosine"
    
    # Minimum chain length (number of carbon atoms) for both the acyl chain and the sphingoid tail.
    MIN_CHAIN_LENGTH = 8

    # Helper function to check that an atom has a hydroxyl group.
    # It now checks not only direct O neighbors (single‐bonded oxygen) but also
    # looks one bond beyond if the neighbor is carbon that is bound to an -OH.
    def has_hydroxyl(atom, core_atoms):
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            # Check for a directly attached oxygen via a single bond.
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    return True
            # If the neighbor is a carbon, check if that carbon has an -OH (i.e. CH2OH).
            if nbr.GetSymbol() == 'C':
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetIdx() == atom.GetIdx() or subnbr.GetIdx() in core_atoms:
                        continue
                    if subnbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            return True
        return False

    # DFS helper to find the longest contiguous carbon chain (only considers single bonds
    # between carbon atoms and avoids atoms in core_atoms).
    def longest_C_chain(atom, visited, core_atoms):
        if atom.GetSymbol() != 'C':
            return 0
        max_len = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited or nbr.GetIdx() in core_atoms:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            if nbr.GetSymbol() != 'C':
                continue
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            chain_len = 1 + longest_C_chain(nbr, new_visited, core_atoms)
            if chain_len > max_len:
                max_len = chain_len
        return max_len

    # Loop over each candidate core match and further verify criteria.
    for match in matches:
        # match is a tuple of 5 atom indices corresponding to:
        #   idx0: acyl carbon (the carbonyl carbon; should have a double-bonded O)
        #   idx1: amide nitrogen
        #   idx2: first head carbon (typically CH2OH)
        #   idx3: second head carbon (bearing -OH)
        #   idx4: third head carbon (bearing -OH and with the sphingoid tail attached)
        acylC = mol.GetAtomWithIdx(match[0])
        amideN = mol.GetAtomWithIdx(match[1])
        head1 = mol.GetAtomWithIdx(match[2])
        head2 = mol.GetAtomWithIdx(match[3])
        head3 = mol.GetAtomWithIdx(match[4])
        
        core_atoms = set(match)
        
        # Check that the acyl carbon has a double-bonded oxygen (the carbonyl).
        has_carbonyl = False
        for nbr in acylC.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(acylC.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue  # try next match
        
        # Verify that each head carbon carries a hydroxyl (either directly or via a CH2OH group).
        if not has_hydroxyl(head1, core_atoms):
            continue
        if not has_hydroxyl(head2, core_atoms):
            continue
        if not has_hydroxyl(head3, core_atoms):
            continue
        
        # Evaluate the fatty acyl chain attached to the acyl carbon.
        acyl_chain_max = 0
        for nbr in acylC.GetNeighbors():
            if nbr.GetIdx() == amideN.GetIdx():
                continue
            visited = {acylC.GetIdx()}
            chain_len = longest_C_chain(nbr, visited, core_atoms)
            if chain_len > acyl_chain_max:
                acyl_chain_max = chain_len
        if acyl_chain_max < MIN_CHAIN_LENGTH:
            continue  # acyl chain too short
        
        # Evaluate the sphingoid tail attached to the third head carbon.
        tail_chain_max = 0
        for nbr in head3.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            visited = {head3.GetIdx()}
            chain_len = longest_C_chain(nbr, visited, core_atoms)
            if chain_len > tail_chain_max:
                tail_chain_max = chain_len
        if tail_chain_max < MIN_CHAIN_LENGTH:
            continue  # sphingoid tail too short
        
        # Passed all criteria for this candidate.
        return True, "Molecule contains an acylated phytosphingosine core with sufficiently long fatty acyl and sphingoid chains"
    
    # If no candidate passed all tests.
    return False, "Matches found for the core motif but did not meet acyl or sphingoid chain criteria"

# For informal testing, you might run:
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
# print(is_N_acylphytosphingosine(test_smiles))
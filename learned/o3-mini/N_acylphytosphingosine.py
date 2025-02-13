"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
The heuristic looks for a core substructure that is roughly:
     –R–C(=O)–N–C(CO)–C(O)–C(O)–
with the expectation that the acyl portion (to the left of the carbonyl)
and the sphingoid chain (attached to the terminal head carbon) each are long aliphatic chains.
Extra checks ensure that each “head” carbon carries a hydroxyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    
    A N-acylphytosphingosine is (roughly) defined as:
         R–C(=O)–N–C(CO)–C(O)–C(O)–[tail]
    where the acyl (fatty acyl) chain (R–) and the sphingoid tail (attached to the third
    head carbon) are long (here at least 8 carbon atoms long) and the three carbon head bears hydroxyl groups.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if molecule is N-acylphytosphingosine, False otherwise.
       str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core SMARTS.
    # It looks for an acyl carbonyl [C](=O) attached to an amide N, which in turn is bound to three carbons:
    #   first carbon with a -CO group (typically CH2OH), then two carbons (each with an -OH).
    core_smarts = "[C](=O)[N][C](CO)[C](O)[C](O)"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return None, None  # Should not happen
    
    matches = mol.GetSubstructMatches(core_pattern)
    if not matches:
        return False, f"Core motif {core_smarts} not found: not an N-acylphytosphingosine"
    
    # Helper: check whether an atom has an -OH substituent (an oxygen attached by a single bond)
    def has_hydroxyl(atom, core_atoms):
        for nbr in atom.GetNeighbors():
            # Only count noncore substituents
            if nbr.GetIdx() in core_atoms:
                continue
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    return True
        return False

    # Helper DFS that returns the longest contiguous chain length (number of carbon atoms)
    # Only traverses atoms that are carbon, not in core_atoms, and only through single bonds.
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

    # We require the acyl chain and sphingoid tail to have at least this many carbon atoms.
    MIN_CHAIN_LENGTH = 8

    # For each candidate core match check the further criteria:
    for match in matches:
        # match is a tuple of 5 atom indices corresponding to:
        #   idx0: acyl carbon (carbonyl carbon)
        #   idx1: amide nitrogen
        #   idx2: first head carbon (should be CH2OH)
        #   idx3: second head carbon (with OH)
        #   idx4: third head carbon (with OH; tail attaches here)
        acylC = mol.GetAtomWithIdx(match[0])
        amideN = mol.GetAtomWithIdx(match[1])
        head1 = mol.GetAtomWithIdx(match[2])
        head2 = mol.GetAtomWithIdx(match[3])
        head3 = mol.GetAtomWithIdx(match[4])
        
        core_atoms = set(match)
        
        # Verify that the acyl carbon (idx0) has a carbonyl oxygen (O=)
        has_carbonyl = False
        for nbr in acylC.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(acylC.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue  # try next match
        
        # Check that the head atoms have the expected hydroxyl groups.
        if not has_hydroxyl(head1, core_atoms):
            continue
        if not has_hydroxyl(head2, core_atoms):
            continue
        if not has_hydroxyl(head3, core_atoms):
            continue
        
        # Evaluate the fatty acyl chain.
        # From acylC, follow all neighbors except the amide nitrogen.
        acyl_chain_max = 0
        for nbr in acylC.GetNeighbors():
            if nbr.GetIdx() == amideN.GetIdx():
                continue
            visited = {acylC.GetIdx()}
            chain_len = longest_C_chain(nbr, visited, core_atoms)
            if chain_len > acyl_chain_max:
                acyl_chain_max = chain_len
        if acyl_chain_max < MIN_CHAIN_LENGTH:
            continue  # too short acyl chain

        # Evaluate the sphingoid tail.
        # From head3 (the third head carbon) follow neighbors that are not part of the core.
        tail_chain_max = 0
        for nbr in head3.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            visited = {head3.GetIdx()}
            chain_len = longest_C_chain(nbr, visited, core_atoms)
            if chain_len > tail_chain_max:
                tail_chain_max = chain_len
        if tail_chain_max < MIN_CHAIN_LENGTH:
            continue  # too short sphingoid tail
        
        # If we pass all tests for this candidate, we classify as an N-acylphytosphingosine.
        return True, "Molecule contains an acylated phytosphingosine core with sufficiently long fatty acyl and sphingoid chains"
    
    # If no candidate passed, choose an appropriate message.
    return False, "Matches found for the core motif but did not meet acyl or sphingoid chain criteria"

# For informal testing, you can uncomment the lines below:
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
# print(is_N_acylphytosphingosine(test_smiles))
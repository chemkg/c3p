"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
This improved version first finds a core motif roughly corresponding to:
    R–C(=O)–N–C(CO)[C](O)[C](O)
Then it verifies that:
  1. The acyl carbon (attached to the C=O) indeed has a double‐bonded oxygen.
  2. Each of the three “head” carbons (the ones from the sphingosine backbone) bears a hydroxyl group,
     checking both direct –OH attachment or an attached CH2OH group.
  3. Both the fatty acyl chain (attached at the acyl carbon) and the sphingoid tail
     (coming off the third head carbon) are sufficiently long (here, at least 8 carbon atoms),
     where chain searches are now forced to ignore the core motif.
If any of these criteria is not met then the molecule is rejected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    
    A N-acylphytosphingosine is (roughly) defined as a ceramide in which a fatty acyl group
    is attached to a phytosphingosine base having the following core motif:
         R–C(=O)–N–C(CO)[C](O)[C](O)
    In addition, both the fatty acyl chain (R–) and the sphingoid tail (attached to the third head carbon)
    are required to be long (here at least 8 carbon atoms long), and each head carbon must carry a hydroxyl group.

    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if molecule is classified as an N-acylphytosphingosine, False otherwise.
       str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core SMARTS.
    # It captures: acyl carbon (with carbonyl), amide N, then three carbons:
    # the first bearing a –CH2OH, then two secondary carbons with –OH.
    # This SMARTS intentionally ignores chirality.
    core_smarts = "[C](=O)[N][C](CO)[C](O)[C](O)"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return None, None  # Should not happen
    
    matches = mol.GetSubstructMatches(core_pattern)
    if not matches:
        return False, f"Core motif {core_smarts} not found: not an N-acylphytosphingosine"
    
    MIN_CHAIN_LENGTH = 8  # minimum number of contiguous carbon atoms required

    # Helper: check that a given atom has a hydroxyl group attached (either directly or via a CH2OH branch)
    def has_hydroxyl(atom, core_atoms):
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            # Check for directly attached oxygen with a single bond (–OH)
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    return True
            # If neighbor is a carbon, check if that carbon carries an –OH (e.g. CH2OH)
            if nbr.GetSymbol() == 'C':
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetIdx() == atom.GetIdx() or subnbr.GetIdx() in core_atoms:
                        continue
                    if subnbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            return True
        return False

    # Helper: Perform a DFS to find the longest contiguous carbon chain,
    # only following single bonds between carbon atoms and ignoring any atoms in core_atoms.
    def longest_C_chain(atom, visited):
        if atom.GetSymbol() != 'C':
            return 0
        max_len = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            if nbr.GetSymbol() != 'C':
                continue
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            chain_len = 1 + longest_C_chain(nbr, new_visited)
            if chain_len > max_len:
                max_len = chain_len
        return max_len

    # Loop over each candidate (the molecule might have multiple substructure matches).
    for match in matches:
        # match is a tuple of 5 atom indices corresponding to:
        #   match[0]: acyl carbon (the carbonyl carbon; should have a double-bonded O)
        #   match[1]: amide nitrogen
        #   match[2]: first head carbon (usually CH2OH)
        #   match[3]: second head carbon (bearing –OH)
        #   match[4]: third head carbon (bearing –OH and with the sphingoid tail attached)
        acylC = mol.GetAtomWithIdx(match[0])
        amideN = mol.GetAtomWithIdx(match[1])
        head1 = mol.GetAtomWithIdx(match[2])
        head2 = mol.GetAtomWithIdx(match[3])
        head3 = mol.GetAtomWithIdx(match[4])
        
        core_atoms = set(match)
        
        # Verify that the acyl carbon has a double-bonded oxygen (carbonyl group)
        has_carbonyl = False
        for nbr in acylC.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(acylC.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue
        
        # Check that each head carbon carries a hydroxyl (either directly or via a CH2OH group)
        if not has_hydroxyl(head1, core_atoms):
            continue
        if not has_hydroxyl(head2, core_atoms):
            continue
        if not has_hydroxyl(head3, core_atoms):
            continue
        
        # Evaluate the fatty acyl chain attached to the acyl carbon.
        acyl_chain_max = 0
        # Among neighbors of the acyl carbon, ignore the one that is the amide nitrogen.
        for nbr in acylC.GetNeighbors():
            if nbr.GetIdx() == amideN.GetIdx():
                continue
            # Start visited as all core atoms so DFS stays out of the core.
            visited = set(core_atoms)
            # Also add the current atom to the visited set.
            visited.add(nbr.GetIdx())
            chain_len = 1 + longest_C_chain(nbr, visited)
            if chain_len > acyl_chain_max:
                acyl_chain_max = chain_len
        if acyl_chain_max < MIN_CHAIN_LENGTH:
            # Explain that the fatty acyl chain is too short.
            continue
        
        # Evaluate the sphingoid tail (should extend from the third head carbon).
        tail_chain_max = 0
        for nbr in head3.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                continue
            visited = set(core_atoms)
            visited.add(nbr.GetIdx())
            chain_len = 1 + longest_C_chain(nbr, visited)
            if chain_len > tail_chain_max:
                tail_chain_max = chain_len
        if tail_chain_max < MIN_CHAIN_LENGTH:
            continue
        
        # All criteria passed for this candidate.
        return True, ("Molecule contains an acylated phytosphingosine core with a carbonyl, "
                      "appropriately hydroxylated head groups, and sufficiently long fatty acyl and sphingoid chains")
    
    # If no candidate passed all the tests.
    return False, "Matches found for the core motif but did not meet acyl or sphingoid chain criteria"

# For example, for testing you might run:
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
# print(is_N_acylphytosphingosine(test_smiles))
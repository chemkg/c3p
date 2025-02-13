"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine

Definition:
An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl groups 
are located at positions 1 and 2 respectively. These molecules consist of three parts:
  • A phosphocholine headgroup (typically contains the substructure COP([O-])(=O)OCC[N+](C)(C)C)
  • One acyl ester chain (an O–C(=O) linkage attached via the glycerol backbone)
  • One O-alkyl (ether) chain that is not part of an ester or phosphate substructure and 
    shows a long uninterrupted aliphatic chain.
    
This program improves on previous attempts by “anchoring” the acyl group to the glycerol 
backbone (by checking that its oxygen is attached to a carbon that neighbors a phosphorus atom)
and by filtering out oxygens that are directly linked to phosphorus when searching for a long O-alkyl chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_alkyl_chain(atom, visited=None):
    """
    Recursively finds the longest chain (in terms of number of carbon-carbon bonds)
    starting from the given carbon atom. Only aliphatic sp3 carbons that are not part
    of a carbonyl (double-bonded to oxygen) are traversed.
    
    Args:
        atom: RDKit Atom from which to start the search.
        visited: Set of atom indices already visited.
    
    Returns:
        int: Length (number of bonds) of the longest alkyl chain found.
    """
    if visited is None:
        visited = set()
    visited.add(atom.GetIdx())
    max_length = 0
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
            # Skip atoms that are part of a carbonyl group.
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
    
    The heuristic consists of:
      1. Checking for the phosphocholine headgroup (using the SMARTS "COP([O-])(=O)OCC[N+](C)(C)C").
      2. Looking for a single acyl ester group [O–C(=O)] that appears to be attached to a glycerol backbone.
         We filter these by requiring that the oxygen in the ester group is attached (through a carbon)
         to an atom that in turn is bonded to phosphorus.
      3. Identifying an O-alkyl (ether) linkage by skipping oxygens that are part of the phosphate
         headgroup (or already used in the ester) and confirming that at least one of its carbon neighbors
         leads to a long uninterrupted alkyl chain (chain length >= 8).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine.
        str: A reason string explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the phosphocholine headgroup.
    # Using a SMARTS pattern that captures the typical phosphocholine motif.
    phos_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phos_pattern):
        return False, "Missing phosphocholine headgroup (pattern COP([O-])(=O)OCC[N+](C)(C)C not found)"
    
    # 2. Look for the acyl ester group pattern.
    # We use the SMARTS "O[C](=O)" to capture an ester oxygen.
    ester_pattern = Chem.MolFromSmarts("O[C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "Missing acyl ester group (pattern O[C](=O) not found)"
    
    # Filter ester matches so that we only count those whose oxygen appears linked to glycerol,
    # i.e. its non-carbonyl neighbor leads (within one bond) to a phosphorus atom.
    valid_ester_count = 0
    ester_oxygen_indices = set()
    for match in ester_matches:
        # In the match, match[0] is the oxygen.
        oxy = mol.GetAtomWithIdx(match[0])
        for nbr in oxy.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # candidate for the glycerol carbon
                # Check if any neighbor of this carbon is phosphorus.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 15:
                        valid_ester_count += 1
                        ester_oxygen_indices.add(oxy.GetIdx())
                        break
                if oxy.GetIdx() in ester_oxygen_indices:
                    break  # found a phosphorus connection; no need to check further
    if valid_ester_count != 1:
        return False, f"Expected exactly one acyl ester group anchored to glycerol but found {valid_ester_count}."
    
    # 3. Identify an O-alkyl (ether) chain.
    # We look for an oxygen that is not part of the ester group and not directly connected to phosphorus.
    ether_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # consider oxygen atoms only
        # Skip oxygens that are already assigned to the ester group.
        if atom.GetIdx() in ester_oxygen_indices:
            continue
        # Skip oxygens that are directly attached to phosphorus.
        if any(neigh.GetAtomicNum() == 15 for neigh in atom.GetNeighbors()):
            continue
        
        # For each remaining oxygen, check among its carbon neighbors for a long alkyl chain.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            chain_length = longest_alkyl_chain(nbr)
            if chain_length >= 8:
                ether_found = True
                break
        if ether_found:
            break
    if not ether_found:
        return False, ("Missing or too short O-alkyl (ether) chain; no oxygen found (outside the ester and headgroup) "
                       "that leads to an alkyl chain of sufficient length (>=8 bonds).")
    
    return True, ("Contains a phosphocholine headgroup, exactly one acyl ester group anchored to glycerol, "
                   "and one O-alkyl (ether) chain with a long uninterrupted alkyl chain. "
                   "This is consistent with 2-acyl-1-alkyl-sn-glycero-3-phosphocholine.")

# For testing purposes (example molecule from provided data):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCC"
    result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(test_smiles)
    print(f"Test result: {result}\nReason: {reason}")
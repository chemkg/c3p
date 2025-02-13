"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: endocannabinoid
A class of cannabinoids present in mammalian biological fluids and tissues that activate cannabinoid receptors.
This improved version uses the following heuristics:
  1. The molecule must contain at least 16 carbon atoms.
  2. It must include either:
       (a) An acyl ethanolamide motif (SMARTS "[C](=O)[N]CCO") in which the acyl chain (attached to the carbonyl)
           is long (≥8 contiguous acyclic carbon atoms), or
       (b) A glycerol backbone motif (using stereospecific SMARTS patterns "O[C@H](CO)CO" or "O[C@@H](CO)CO")
           where none of the glycerol atoms are in rings and at least one glycerol oxygen is attached to a long acyclic
           aliphatic chain.
Improvements include restrictions so that only acyclic chains (thus avoiding aromatic or ring‐embedded carbons)
are counted and ensuring that the glycerol backbone is not part of a fused ring system.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def _longest_chain(mol, start_idx, visited):
    """
    Recursively computes the length of a contiguous acyclic carbon chain starting from a given atom.
    Only carbon atoms that are not in rings are traversed.
    """
    atom = mol.GetAtomWithIdx(start_idx)
    max_length = 1
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # Only traverse carbon atoms that are not in a ring
        if nbr.GetAtomicNum() == 6 and not nbr.IsInRing() and nbr_idx not in visited:
            new_visited = visited.union({nbr_idx})
            chain_length = 1 + _longest_chain(mol, nbr_idx, new_visited)
            if chain_length > max_length:
                max_length = chain_length
    return max_length

def _check_long_acyl(mol, start_idx, min_length=8):
    """
    Check if starting from the given carbon atom there is a contiguous acyclic chain of at least min_length carbon atoms.
    """
    return _longest_chain(mol, start_idx, {start_idx}) >= min_length

def _glycerol_non_ring(mol, atom_indices):
    """
    Verifies that all atoms in the specified list (a candidate glycerol match) are not part of a ring.
    """
    for idx in atom_indices:
        if mol.GetAtomWithIdx(idx).IsInRing():
            return False
    return True

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string using improved heuristics.
    
    Heuristic criteria:
      1. The molecule must have at least 16 carbon atoms.
      2. It must contain either:
         (a) An acyl ethanolamide motif (SMARTS "[C](=O)[N]CCO")
             in which the acyl chain (starting from the carbonyl carbon’s carbon neighbor) is long (≥8 contiguous acyclic carbons),
             or
         (b) A glycerol backbone motif defined by stereospecific patterns ("O[C@H](CO)CO" or "O[C@@H](CO)CO"),
             with the extra requirement that the glycerol match is entirely non‐cyclic.
             In addition, at least one glycerol oxygen must be linked to a long acyclic aliphatic chain.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule meets the endocannabinoid heuristic, False otherwise.
      str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Filter 1: Overall carbon count must be at least 16.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, "Too few carbon atoms to support a long fatty acyl chain typical of endocannabinoids"
    
    # (a) Check for acyl ethanolamide motif using SMARTS "[C](=O)[N]CCO"
    ethanolamide_pattern = Chem.MolFromSmarts("[C](=O)[N]CCO")
    if mol.HasSubstructMatch(ethanolamide_pattern):
        matches = mol.GetSubstructMatches(ethanolamide_pattern)
        for match in matches:
            # In the match, the first atom is the carbonyl carbon.
            carbonyl_idx = match[0]
            carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
            # Look for a carbon neighbor (excluding the amide nitrogen) that could be the start of the acyl chain.
            for nbr in carbonyl_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                    if _check_long_acyl(mol, nbr.GetIdx(), min_length=8):
                        return True, ("Contains an acyl ethanolamide moiety with a long fatty acyl chain, "
                                      "characteristic of many endocannabinoids.")
    
    # (b) Check for a glycerol backbone using stereospecific SMARTS patterns.
    glycerol_smarts_list = ["O[C@H](CO)CO", "O[C@@H](CO)CO"]
    for gs in glycerol_smarts_list:
        glycerol_pattern = Chem.MolFromSmarts(gs)
        if mol.HasSubstructMatch(glycerol_pattern):
            matches = mol.GetSubstructMatches(glycerol_pattern)
            for match in matches:
                # Exclude matches where any glycerol atom is part of a ring.
                if not _glycerol_non_ring(mol, match):
                    continue
                # Now check: at least one oxygen in the glycerol match should be linked to a long acyclic carbon chain.
                for atom_idx in match:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == "O":
                        for nbr in atom.GetNeighbors():
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in match and not nbr.IsInRing():
                                if _check_long_acyl(mol, nbr.GetIdx(), min_length=8):
                                    return True, ("Contains a glycerol backbone (non-cyclic) with an attached long acyclic "
                                                  "aliphatic chain (via ester or ether linkage), typical of monoacylglycerol endocannabinoids.")
    
    return False, "Does not contain characteristic endocannabinoid structural motifs."

# Example usage (for testing purposes):
# test_smiles = "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)NCCO"  # anandamide example
# result, explanation = is_endocannabinoid(test_smiles)
# print(result, explanation)
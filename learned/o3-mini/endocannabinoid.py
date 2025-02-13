"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: endocannabinoid
A class of cannabinoids present in mammalian biological fluids and tissues that activate cannabinoid receptors.
Improved heuristics:
  1. The molecule must contain at least 16 carbon atoms.
  2. It must contain either:
       (a) An acyl ethanolamide motif – defined by a carbonyl (not in a ring) attached to an amide N and then “CCO.”
           In addition the acyl chain (attached to the carbonyl carbon) must have at least 8 contiguous acyclic carbon atoms.
       (b) A monoacylglycerol motif – a stereospecific glycerol backbone match (using patterns "O[C@H](CO)CO" or "O[C@@H](CO)CO")
           where none of the glycerol atoms belong to a ring and exactly one of its hydroxyl oxygens is acylated (via an ester linkage)
           with a long acyclic chain (≥8 contiguous carbons).
If neither motif is found the molecule is not classified as an endocannabinoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def _longest_chain(mol, start_idx, visited):
    """
    Recursively computes the length of a contiguous acyclic carbon chain starting from the given atom.
    Only traverses carbon atoms not in rings.
    """
    atom = mol.GetAtomWithIdx(start_idx)
    max_length = 1
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr.GetAtomicNum() == 6 and (nbr_idx not in visited) and (not nbr.IsInRing()):
            new_visited = visited.union({nbr_idx})
            chain_length = 1 + _longest_chain(mol, nbr_idx, new_visited)
            if chain_length > max_length:
                max_length = chain_length
    return max_length

def _check_long_acyl(mol, start_idx, min_length=8):
    """
    Check if starting from the given carbon atom there is a contiguous acyclic chain of at least min_length carbons.
    """
    return _longest_chain(mol, start_idx, {start_idx}) >= min_length

def _glycerol_non_ring(mol, atom_indices):
    """
    Verifies that none of the atoms in the candidate glycerol match (provided as a list of indices) are part of a ring.
    """
    for idx in atom_indices:
        if mol.GetAtomWithIdx(idx).IsInRing():
            return False
    return True

def _count_glycerol_acylated_oxygens(mol, match):
    """
    For a glycerol match (list of atom indices), count how many oxygen atoms are acylated.
    An oxygen is considered acylated if it is directly bound to a carbon that is double-bonded to an oxygen (i.e. in an ester).
    """
    count = 0
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == "O":
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    # Check if this neighbor carbon has a double bond to a non-ring oxygen (i.e. the carbonyl)
                    for bond in nbr.GetBonds():
                        # We require a double bond and that the other atom is oxygen.
                        if bond.GetBondTypeAsDouble() == 2:
                            other = bond.GetOtherAtom(nbr)
                            if other.GetSymbol() == "O" and (not other.IsInRing()):
                                count += 1
                                break
                    # If found one acylation, no need to count multiple (but could be more if accidentally duplicated)
    return count

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string using improved heuristics.
    
    Heuristic criteria:
      1. The molecule must have at least 16 carbon atoms.
      2. It must contain either:
         (a) An acyl ethanolamide motif – a C(=O)NCCO fragment (with the carbonyl not in a ring)
             where the carbonyl carbon is attached to a long acyclic aliphatic chain (≥8 contiguous C's), or
         (b) A monoacylglycerol motif defined by a stereospecific glycerol pattern ("O[C@H](CO)CO" or "O[C@@H](CO)CO")
             in which all glycerol atoms are non-cyclic and exactly one of its hydroxyl oxygens is acylated with a long acyclic chain.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule meets the endocannabinoid heuristic, False otherwise.
      str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Requirement 1: At least 16 carbons overall.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, "Too few carbon atoms to support a long fatty acyl chain typical of endocannabinoids"
    
    # Heuristic (a): Acyl ethanolamide motif.
    # Use a relaxed SMARTS pattern for the motif – note that the carbonyl should not be in a ring.
    ethanolamide_pattern = Chem.MolFromSmarts("[C;!R](=O)[N]CCO")
    if ethanolamide_pattern is not None and mol.HasSubstructMatch(ethanolamide_pattern):
        matches = mol.GetSubstructMatches(ethanolamide_pattern)
        for match in matches:
            # match[0]: carbonyl carbon
            # match[1]: amide nitrogen; match[2]: CH2 next to nitrogen; match[3]: CH2-O part of ethanolamine.
            carbonyl_atom = mol.GetAtomWithIdx(match[0])
            if carbonyl_atom.IsInRing():
                continue  # skip if in a ring
            # Find carbon neighbor(s) of the carbonyl that are not the amide nitrogen.
            for nbr in carbonyl_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in match:
                    if not nbr.IsInRing() and _check_long_acyl(mol, nbr.GetIdx(), min_length=8):
                        return True, ("Contains an acyl ethanolamide moiety with a long fatty acyl chain, "
                                      "characteristic of many endocannabinoids.")
    
    # Heuristic (b): Monoacylglycerol motif.
    # Look for glycerol backbone stereospecific patterns.
    glycerol_smarts_list = ["O[C@H](CO)CO", "O[C@@H](CO)CO"]
    for gs in glycerol_smarts_list:
        glycerol_pattern = Chem.MolFromSmarts(gs)
        if glycerol_pattern is None:
            continue
        if mol.HasSubstructMatch(glycerol_pattern):
            matches = mol.GetSubstructMatches(glycerol_pattern)
            for match in matches:
                # Ensure none of the glycerol atoms in this match are in rings.
                if not _glycerol_non_ring(mol, match):
                    continue
                # Count acylated oxygens on this glycerol fragment.
                acylated_count = _count_glycerol_acylated_oxygens(mol, match)
                # For a monoacylglycerol endocannabinoid, exactly one hydroxyl (oxygen) should be acylated.
                if acylated_count != 1:
                    continue
                # Identify the acylated oxygen and then check the attached carbonyl’s neighbor for a long chain.
                for idx in match:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() != "O":
                        continue
                    for nbr in atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 6:
                            # Look for a bond where nbr acts as the carbonyl carbon (has a double bond to oxygen)
                            for bond in nbr.GetBonds():
                                if bond.GetBondTypeAsDouble() == 2:
                                    other = bond.GetOtherAtom(nbr)
                                    if other.GetSymbol() == "O":
                                        # Found an ester linkage.
                                        # Now find the neighbor of the carbonyl (nbr) that is not the ester oxygen.
                                        for nb in nbr.GetNeighbors():
                                            if nb.GetIdx() == idx:
                                                continue
                                            if nb.GetAtomicNum() == 6 and not nb.IsInRing():
                                                if _check_long_acyl(mol, nb.GetIdx(), min_length=8):
                                                    return True, ("Contains a monoacylglycerol motif: a non‐cyclic glycerol backbone acylated "
                                                                  "at one hydroxyl with a long aliphatic chain, typical of endogenous cannabinoids.")
    return False, "Does not contain characteristic endocannabinoid structural motifs."

# Example usage (for testing purposes):
# test_smiles = "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)NCCO"  # anandamide example
# result, explanation = is_endocannabinoid(test_smiles)
# print(result, explanation)
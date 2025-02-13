"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: endocannabinoid
A class of cannabinoids present in mammalian biological fluids and tissues that activate cannabinoid receptors.
This program uses improved structural heuristics:
  1. The molecule must have enough carbons (here ≥ 16) to support a fatty chain.
  2. The molecule must include at least one of two motifs:
       (a) An acyl ethanolamide motif (e.g. anandamide, palmitoyl ethanolamide) in which the acyl group is long (≥8 carbons).
       (b) A glycerol backbone (using stereochemically defined SMARTS) that has an attached long aliphatic chain (via an ester or ether linkage),
           found in monoacylglycerol endocannabinoids (and extended here to cover both ester and ether linkages).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def _longest_chain(mol, start_idx, visited):
    """
    Recursively determines the length of the contiguous carbon chain starting at a given atom
    by only traversing carbon atoms.
    """
    atom = mol.GetAtomWithIdx(start_idx)
    max_length = 1
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:  # only carbons
            new_visited = visited.union({nbr_idx})
            chain_length = 1 + _longest_chain(mol, nbr_idx, new_visited)
            if chain_length > max_length:
                max_length = chain_length
    return max_length

def _check_long_acyl(mol, start_idx, min_length=8):
    """
    Check if starting from the given index (a carbon in the acyl chain)
    one can follow a contiguous path of carbon atoms of at least min_length.
    """
    return _longest_chain(mol, start_idx, {start_idx}) >= min_length

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    
    Improved heuristic criteria:
      1. The molecule must have at least 16 carbon atoms.
      2. It must contain either:
           (a) An acyl ethanolamide motif (SMARTS "[C](=O)[N]CCO") where the carbonyl (acyl) carbon
               leads into a long (≥8 carbons) fatty chain; or
           (b) A glycerol backbone motif – here defined using the stereochemical SMARTS "O[C@H](CO)CO"
               (or its enantiomer) – and at least one of the glycerol oxygens must be linked to a long
               aliphatic chain (the chain may be attached via ester OR ether linkage).
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule meets the endocannabinoid heuristic, False otherwise.
      str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (1) Overall carbon count filter.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, "Too few carbon atoms to support a long fatty acyl chain typical of endocannabinoids"
    
    # (2a) Check for the acyl ethanolamide motif.
    # SMARTS: a carbonyl directly attached to an N followed by –CCO.
    ethanolamide_pattern = Chem.MolFromSmarts("[C](=O)[N]CCO")
    if mol.HasSubstructMatch(ethanolamide_pattern):
        # For each match, verify that the fatty (acyl) chain emerging from the carbonyl is long.
        # In the pattern, the first atom (index 0) is the carbonyl carbon.
        matches = mol.GetSubstructMatches(ethanolamide_pattern)
        for match in matches:
            carbonyl_idx = match[0]
            carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
            # Look at neighbors of the carbonyl carbon that are carbon and not part of the amide linkage.
            for nbr in carbonyl_atom.GetNeighbors():
                # (Skip the neighbor that is the amide nitrogen; typical match order: [C](=O)[N]...)
                if nbr.GetAtomicNum() == 6:
                    if _check_long_acyl(mol, nbr.GetIdx()):
                        return True, ("Contains an acyl ethanolamide moiety with a long fatty acyl chain, "
                                      "characteristic of many endocannabinoids.")
    
    # (2b) Check for a glycerol-based motif.
    # Use a more stereospecific glycerol pattern to reduce false positives.
    glycerol_smarts_list = ["O[C@H](CO)CO", "O[C@@H](CO)CO"]
    for gs in glycerol_smarts_list:
        glycerol_pattern = Chem.MolFromSmarts(gs)
        if mol.HasSubstructMatch(glycerol_pattern):
            matches = mol.GetSubstructMatches(glycerol_pattern)
            for match in matches:
                # The matched atoms are part of a glycerol backbone.
                # Now check each oxygen in the match to see if it is tethered to a long aliphatic chain.
                for atom_idx in match:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == "O":
                        for nbr in atom.GetNeighbors():
                            # Only consider those neighbor carbons that are NOT part of the glycerol match.
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in match:
                                if _check_long_acyl(mol, nbr.GetIdx()):
                                    return True, ("Contains a glycerol backbone with an attached long acyl chain "
                                                  "(via ester or ether linkage), typical of monoacylglycerol endocannabinoids.")
    
    return False, "Does not contain characteristic endocannabinoid structural motifs."

# Example usage (uncomment for testing):
# test_smiles = "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)NCCO"  # anandamide example
# result, reason = is_endocannabinoid(test_smiles)
# print(result, reason)
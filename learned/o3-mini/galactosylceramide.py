"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.
This version improves upon earlier attempts by:
  1. Requiring that exactly one galactose ring (alpha or beta) is detected.
  2. Scanning through each amide bond to find a candidate acyl chain that has at least 10 consecutive “fatty” carbon atoms,
     where we allow any sp3 carbon that is nonaromatic (i.e. not strictly CH2 but also CH(OH) is permitted) as part of a linear chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType, HybridizationType

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is defined as a glycosphingolipid consisting of a ceramide 
    (a sphingoid base linked via an amide bond to a fatty acid) and a single galactose head group.
    
    This function applies the following filters:
      1. The molecule must contain exactly one galactose ring (either alpha or beta form).
      2. The molecule must contain an amide bond. For at least one such bond the acyl portion must have a candidate
         “linear” fatty acyl chain with at least 10 consecutive aliphatic (sp3, non‐aromatic) carbon atoms.
      3. Crude overall filters on molecular weight and total carbon count are applied.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a galactosylceramide, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Sugar head group check ---
    # Define substructure patterns for beta- and alpha-galactose.
    beta_gal = Chem.MolFromSmiles("CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    alpha_gal = Chem.MolFromSmiles("CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    
    gal_matches = []
    if beta_gal:
        gal_matches.extend(mol.GetSubstructMatches(beta_gal))
    if alpha_gal:
        gal_matches.extend(mol.GetSubstructMatches(alpha_gal))
    # Require exactly one galactose ring match.
    if len(gal_matches) != 1:
        return False, f"Expected exactly one galactose sugar head group, found {len(gal_matches)}"
    
    # --- Fatty acyl chain detection from an amide bond ---
    # We search for a carbonyl carbon (C=O) that is part of an amide.
    fatty_acyl_found = False
    
    # Define a helper that returns True if an atom is an acceptable candidate for being part of an acyl chain.
    def is_fatty_chain_atom(atom):
        # Accept only carbon (atomic number 6) that is sp3 hybridized and nonaromatic.
        if atom.GetAtomicNum() != 6:
            return False
        if atom.GetIsAromatic():
            return False
        if atom.GetHybridization() != HybridizationType.SP3:
            return False
        return True

    # Depth-first search to compute the length (number of atoms) in a linear fatty chain.
    def dfs_chain_length(atom, prev, visited):
        max_length = 0
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == prev.GetIdx():
                continue
            # Only continue if the neighbor is an acceptable fatty chain atom.
            if is_fatty_chain_atom(nbr) and (nbr.GetIdx() not in visited):
                length = 1 + dfs_chain_length(nbr, atom, visited)
                if length > max_length:
                    max_length = length
        # Remove the atom from visited (to allow other paths in different branches)
        visited.remove(atom.GetIdx())
        return max_length

    # Loop over all atoms and look for a carbonyl carbon that is part of an amide.
    for atom in mol.GetAtoms():
        # Look for a carbon atom that is double-bonded to oxygen.
        if atom.GetAtomicNum() != 6:
            continue
        carbonyl_found = False
        amide_nitrogen = False
        oxygen_neighbor = None
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            # Check for a double bond to oxygen.
            if bond.GetBondType() == BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                carbonyl_found = True
                oxygen_neighbor = nbr
            # Check for a single bond to a nitrogen.
            if bond.GetBondType() == BondType.SINGLE and nbr.GetAtomicNum() == 7:
                amide_nitrogen = True
        if carbonyl_found and amide_nitrogen:
            # This atom is a candidate for an amide carbonyl.
            # Identify the neighbor which is not the carbonyl oxygen and not the amide nitrogen,
            # which should lead to the fatty acyl chain.
            acyl_neighbor = None
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == (oxygen_neighbor.GetIdx() if oxygen_neighbor else -1):
                    continue
                if nbr.GetAtomicNum() == 7:  # Skip nitrogen
                    continue
                acyl_neighbor = nbr
                break
            # If an acyl chain candidate exists and is a fatty chain atom (or even if not perfectly CH2),
            # do a DFS to find the longest linear chain.
            if acyl_neighbor and is_fatty_chain_atom(acyl_neighbor):
                chain_len = 1 + dfs_chain_length(acyl_neighbor, atom, set())
                if chain_len >= 10:
                    fatty_acyl_found = True
                    break
    if not fatty_acyl_found:
        return False, "No long fatty acyl chain (>=10 consecutive acceptable aliphatic carbons) detected on an amide bond – not a typical ceramide"
    
    # --- Crude overall filters ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a typical glycosphingolipid or ceramide derivative"
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms to represent a typical ceramide derivative"
    
    return True, ("Contains exactly one galactose sugar head group and an amide-linked ceramide "
                  "backbone with a long fatty acyl chain (galactosylceramide)")

# Example tests (uncomment to run):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_galactosylceramide(test_smiles)
# print(result, reason)
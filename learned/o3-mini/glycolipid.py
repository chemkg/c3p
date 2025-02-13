"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid

Definition: A glycolipid is defined as a molecule that comprises a carbohydrate (sugar)
moiety linked glycosidically to a lipid part. Typical examples include 1,2-di-O-acylglycerols 
linked to one mono-, di-, or trisaccharide, as well as glycosphingolipids and bacterial glycolipids 
that may feature acylated sugar parts and/or an amide-linked long acyl chain.

Heuristic improvements made:
 • The sugar-detection helper now requires a 5- or 6-membered ring containing at least one oxygen 
   and at least one exocyclic hydroxyl group (rather than two) so that partially acylated sugars are not missed.
 • In addition, if no sugar ring is found, we also try a simple substructure search for a putative 
   anomeric carbon pattern (to catch cases where acylation may block many –OH groups).
 • The lipid-detection heuristics remain based on searching for long-chain ester/amide groups or a long 
   aliphatic chain attached to the sugar moiety.
 
This remains an approximate classification.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_sugar_rings(mol):
    """
    Identify sugar-like rings as 5- or 6-membered rings that contain at least one oxygen
    and at least one exocyclic hydroxyl group (an oxygen outside the ring bonded to hydrogen).
    Returns a list of sets, each a set of atom indices.
    """
    sugar_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Must contain at least one oxygen in the ring.
            if not any(atom.GetSymbol() == "O" for atom in ring_atoms):
                continue
            hydroxyl_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    # Exclude neighbors that are also in the ring.
                    if nbr.GetIdx() in ring:
                        continue
                    # Look for exocyclic oxygen that is bonded to at least one hydrogen.
                    if nbr.GetSymbol() == "O":
                        if any(nbr2.GetSymbol() == "H" for nbr2 in nbr.GetNeighbors()):
                            hydroxyl_count += 1
            # Relaxed threshold: require at least one exocyclic hydroxyl instead of two.
            if hydroxyl_count >= 1:
                sugar_rings.append(set(ring))
    return sugar_rings

def search_for_anomeric_pattern(mol):
    """
    As a fallback, search for an anomeric carbon pattern often found in sugars.
    Typical pattern: a tetrahedral carbon bonded to two oxygens (one exocyclic and one in a ring).
    This is not perfect but may help detect sugar-like moieties when the ring hydroxyl count is low.
    """
    # This SMARTS looks for a tetrahedral carbon (sp3, nonaromatic) that is bonded to an oxygen
    # that is itself bonded to a ring oxygen.
    anomeric_smarts = "[C;!R;!$(C=O)]-O-[O;R]"
    pattern = Chem.MolFromSmarts(anomeric_smarts)
    return mol.HasSubstructMatch(pattern)

def dfs_long_chain(mol, current_idx, visited, excluded):
    """
    Recursively determines the longest contiguous chain (number of carbon atoms)
    starting from an atom (current_idx). Only considers aliphatic (sp3 nonaromatic) carbons.
    'visited' avoids cycles; 'excluded' is a set of atoms (for example sugar atoms) we skip.
    Returns the chain length counting the starting atom as 1.
    """
    current_atom = mol.GetAtomWithIdx(current_idx)
    max_len = 1
    for nbr in current_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx in visited or nbr_idx in excluded:
            continue
        if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() and nbr.GetHybridization().name == "SP3":
            visited.add(nbr_idx)
            chain_len = 1 + dfs_long_chain(mol, nbr_idx, visited, excluded)
            visited.remove(nbr_idx)
            if chain_len > max_len:
                max_len = chain_len
    return max_len

def longest_aliphatic_chain_from(mol, start_idx, excluded):
    """
    Computes the maximum chain length starting from start_idx over contiguous aliphatic carbons.
    """
    return dfs_long_chain(mol, start_idx, {start_idx}, excluded)

def has_attached_long_chain(mol, sugar_ring):
    """
    For all oxygen atoms directly attached to atoms in the sugar ring (but not part of the ring),
    check whether they are linked (through an O—C bond) to a carbon with a long chain (≥ 6 carbons).
    """
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == "O" and nbr.GetIdx() not in sugar_ring:
                for o_nbr in nbr.GetNeighbors():
                    if o_nbr.GetIdx() == idx:
                        continue
                    if o_nbr.GetAtomicNum() == 6 and not o_nbr.GetIsAromatic() and o_nbr.GetHybridization().name == "SP3":
                        chain_length = longest_aliphatic_chain_from(mol, o_nbr.GetIdx(), excluded=sugar_ring)
                        if chain_length >= 6:
                            return True
    return False

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    
    Heuristic improved criteria:
      • Sugar moiety: Preferably a sugar-like ring (5- or 6-membered with ≥1 oxygen and at least one free -OH)
        OR if not found, a fallback substructure (anomeric pattern) is checked.
      • Lipid part: Presence of a long-chain ester group ([CX3](=O)OCCCCCC) or an amide (NC(=O)CCCCCC)
        or a long contiguous aliphatic chain (≥8 carbons), especially if attached to a sugar ring.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if molecule appears to be a glycolipid, False otherwise.
      str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Identify sugar moieties ---
    sugar_rings = get_sugar_rings(mol)
    if not sugar_rings:
        # fallback: try to detect anomeric carbon pattern common in sugars
        if not search_for_anomeric_pattern(mol):
            return False, "No clear carbohydrate (sugar) moiety detected"
    
    # --- Step 2: Check for lipid features ---
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)OCCCCCC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    n_ester = len(ester_matches)
    
    amide_pattern = Chem.MolFromSmarts("NC(=O)CCCCCC")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    
    overall_long_chain = False
    max_chain = 0
    # Exclude sugar atoms when searching for a long chain.
    sugar_atoms = set()
    for ring in sugar_rings:
        sugar_atoms.update(ring)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and atom.GetHybridization().name == "SP3":
            if atom.GetIdx() in sugar_atoms:
                continue
            chain_len = longest_aliphatic_chain_from(mol, atom.GetIdx(), excluded=sugar_atoms)
            if chain_len > max_chain:
                max_chain = chain_len
    if max_chain >= 8:
        overall_long_chain = True
        
    # Also check for a lipid chain directly attached to a sugar ring.
    attached_chain = False
    for ring in sugar_rings:
        if has_attached_long_chain(mol, ring):
            attached_chain = True
            break

    # --- Step 3: Decision logic ---
    if n_ester >= 2:
        return True, f"Contains sugar moiety and {n_ester} long-chain ester group(s) (diacyl structure expected)"
    elif n_amide >= 1 and overall_long_chain:
        return True, f"Contains sugar moiety and an amide bond with a long acyl chain"
    elif attached_chain:
        return True, "Contains sugar moiety and a sugar-linked long aliphatic chain (glycosidic lipid variant)"
    elif overall_long_chain:
        return True, "Contains sugar moiety and a long aliphatic chain, possible glycolipid variant"
    else:
        return False, "Sugar is present but no clear lipid (ester, amide-linked, or attached long-chain) feature detected"

# Example usage:
if __name__ == '__main__':
    # A simplified glycolipid-like molecule for demonstration.
    test_smiles = "OC1C(O)C(O)C(O)C(O)C1OCC(=O)CCCCCCCC"
    classification, reason = is_glycolipid(test_smiles)
    print("Is glycolipid?", classification)
    print("Reason:", reason)
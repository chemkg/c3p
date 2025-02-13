"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid

Definition: A glycolipid is defined as a molecule that comprises a carbohydrate (sugar)
moiety glycosidically linked to a lipid part. Typical examples include 1,2-di-O-acylglycerols
linked to one mono-, di-, or trisaccharide, as well as glycosphingolipids and bacterial glycolipids
that may feature acylated sugar parts and/or an amide-linked long acyl chain.
 
Heuristic improvements:
 • Use a helper to scan for sugar rings (ring size 5 or 6 with ≥1 oxygen and ≥2 hydroxyl substituents).
 • Detect lipid parts by (1) searching for ester groups [CX3](=O)OCCCCCC and amide groups NC(=O)CCCCCC,
   and (2) for any sugar ring, look at exocyclic oxygen atoms that may connect to a long aliphatic chain.
 • Compute the longest contiguous non-aromatic (sp³) carbon chain from a candidate attachment.
 
This is an approximate classification.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_sugar_rings(mol):
    """
    Returns a list of sets, each a set of atom indices for ring(s) that appear sugar‐like.
    A ring is considered sugar-like if its size is 5 or 6, contains at least one oxygen, and
    at least two substituents that look like hydroxyl (an exocyclic oxygen with at least one hydrogen).
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
                # check neighbors outside of the ring for OH groups
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring: 
                        continue
                    if nbr.GetSymbol() == "O":
                        # check if the oxygen has at least one hydrogen neighbor
                        if any(nbr2.GetSymbol() == "H" for nbr2 in nbr.GetNeighbors()):
                            hydroxyl_count += 1
            if hydroxyl_count >= 2:
                sugar_rings.append(set(ring))
    return sugar_rings

def dfs_long_chain(mol, current_idx, visited, excluded):
    """
    Recursively determine the longest contiguous chain (number of carbon atoms)
    starting from an atom 'current_idx'. Only considers sp3, nonaromatic carbons.
    'visited' is used to avoid cycles.
    'excluded' is a set of atom indices that should not be used (for example, sugar atoms).
    Returns the chain length (counting the starting atom as 1).
    """
    current_atom = mol.GetAtomWithIdx(current_idx)
    max_len = 1
    for nbr in current_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # only continue if neighbor is a carbon, sp3 (nonaromatic) and not in excluded or visited.
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
    check whether they are linked (through an O—C bond) to a carbon whose aliphatic chain has length >= 6.
    """
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # if neighbor is an oxygen and is not in the sugar ring, then check its neighbors:
            if nbr.GetSymbol() == "O" and nbr.GetIdx() not in sugar_ring:
                # look at neighbors of this oxygen (besides our sugar atom)
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
    
    Improved heuristic:
     • Detects a carbohydrate moiety by searching for sugar-like rings.
     • Checks for a lipid part by either:
         - Finding at least two ester groups ([CX3](=O)OCCCCCC) or one amide bond (NC(=O)CCCCCC),
         - OR by locating a sugar exocyclic oxygen that is linked to a long aliphatic chain (chain of ≥6 carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets glycolipid criteria, False otherwise.
        str: Description for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Identify sugar-like rings ---
    sugar_rings = get_sugar_rings(mol)
    if not sugar_rings:
        return False, "No clear carbohydrate (sugar) moiety detected"
    
    # --- Step 2: Check for lipid features ---
    # Search for ester groups with a long acyl chain.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)OCCCCCC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    n_ester = len(ester_matches)
    
    # Search for amide groups with a long acyl chain.
    amide_pattern = Chem.MolFromSmarts("NC(=O)CCCCCC")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    
    # Additionally, search for long aliphatic chains in the overall molecule.
    # Instead of a simple SMILES substring, we compute the longest chain.
    overall_long_chain = False
    max_chain = 0
    for atom in mol.GetAtoms():
        # only start from sp3 carbon atoms not in any sugar ring
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and atom.GetHybridization().name == "SP3":
            # be cautious: only start if not already clearly in a sugar ring
            in_sugar = any(atom.GetIdx() in ring for ring in sugar_rings)
            if in_sugar:
                continue
            chain_len = longest_aliphatic_chain_from(mol, atom.GetIdx(), excluded=set())
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
    # We require that a sugar is present and one of the following lipid features is observed.
    if n_ester >= 2:
        return True, f"Contains sugar moiety and {n_ester} long-chain ester group(s) (diacyl structure expected)"
    elif n_amide >= 1 and overall_long_chain:
        return True, f"Contains sugar moiety and an amide bond with a long acyl chain"
    elif attached_chain:
        return True, "Contains sugar moiety and a sugar-linked long aliphatic chain, possible glycolipid variant"
    elif overall_long_chain:
        # fallback if a very long chain exists
        return True, "Contains sugar moiety and a long aliphatic chain, possible glycolipid variant"
    else:
        return False, "Sugar is present but no clear lipid (ester, amide-linked, or attached long-chain) feature detected"

# Example usage:
if __name__ == '__main__':
    # A simplified glycolipid-like molecule (for demonstration)
    test_smiles = "OC1C(O)C(O)C(O)C(O)C1OCC(=O)CCCCCCCC"
    classification, reason = is_glycolipid(test_smiles)
    print("Is glycolipid?", classification)
    print("Reason:", reason)
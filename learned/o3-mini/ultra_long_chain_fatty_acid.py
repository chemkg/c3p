"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: Ultra Long‐Chain Fatty Acid
Definition: A fatty acid with a free carboxylic acid group (C(=O)[OH]) attached to a nearly linear (or mildly branched)
            carbon chain where the chain (counting the acid carbon) is > 27. In addition, if the molecule contains extra
            substructure (for example large rings outside the fatty acyl fragment) it is less likely to be a “pure” fatty acid.
            
The approach:
  1. Look for a free carboxylic acid group via SMARTS "[CX3](=O)[OX2H1]".
  2. Verify that the acid carbon is attached to exactly one carbon neighbor.
  3. “Walk” from that neighbor to extract one longest linear chain. In doing so we allow walking over carbons that
     are not in rings OR are only found in 3‐membered rings (e.g. cyclopropyl groups).
  4. Form a fatty acid fragment from: (acid carbon + its two oxygens) ∪ (all carbons in the extracted chain).
  5. Also check that there are no “large” rings (rings with >3 atoms) outside the fragment.
  6. If the chain length (number of carbons including the acid carbon) is > 27 and the fatty acid fragment makes up
     “most” of the molecule (here either a purity ratio above about 50% or only very few extra heavy atoms), then
     classify as ultra‐long‐chain fatty acid.
     
Note: This is an approximate algorithm. If the task is too complex to capture all cases, the function may return (None, None).
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines whether a molecule qualifies as an ultra long-chain fatty acid.
    
    The molecule must contain a free carboxylic acid group (C(=O)[OH]) where the carboxyl carbon has exactly one carbon neighbor.
    Starting from that carbon neighbor, we extract a (nearly) linear chain (allowing cyclopropyl groups) and then add the acid group.
    The total carbon chain length (including the acid carbon) must exceed 27.
    In addition, we require that the fatty acid fragment is “pure” in that either the fraction of heavy atoms belonging to the fragment
    is at least 50% or the extra heavy atoms not in the fragment are very few. Finally, if the molecule has any rings bigger than 3 atoms
    that lie outside the fatty acid fragment, we do not classify it as a simple fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an ultra long-chain fatty acid, False otherwise.
        str: Explanation (or reason for rejection).
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    
    # Precompute ring information: list of rings (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Helper: returns True if an atom (by index) belongs to any ring of size >3.
    def in_large_ring(atom_idx):
        for ring in ring_info:
            if len(ring) > 3 and atom_idx in ring:
                return True
        return False
    
    # Helper: Check if an atom is an allowed carbon.
    # Allowed if it is carbon and either not in any ring or only in a cyclopropyl (3-membered) ring.
    def allowed_carbon(atom):
        if atom.GetAtomicNum() != 6:
            return False
        # if atom is in any ring, check if any such ring has size >3:
        idx = atom.GetIdx()
        for ring in ring_info:
            if idx in ring and len(ring) > 3:
                return False
        return True
    
    # Look for free carboxylic acid group: C(=O)[OX2H1]
    acid_smarts = "[CX3](=O)[OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free carboxylic acid group found"
    
    # We consider each acid group separately; if one qualifies then we classify as ultra-long chain.
    reasons = []
    for match in acid_matches:
        # In the SMARTS "[CX3](=O)[OX2H1]", assume index 0 is the acid carbon and index 1 is the hydroxyl oxygen.
        acid_c = mol.GetAtomWithIdx(match[0])
        # The acid carbon should have exactly one carbon neighbor (the fatty acyl chain)
        carbon_neighbors = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            reasons.append(f"Acid carbon (idx {acid_c.GetIdx()}) does not have exactly one carbon neighbor (found {len(carbon_neighbors)})")
            continue
        chain_start = carbon_neighbors[0]
        
        # Recursive DFS to find the longest linear (nonbranching) path
        def get_longest_chain(curr, parent):
            # Only allow if current atom is allowed (using allowed_carbon)
            if not allowed_carbon(curr):
                return [curr.GetIdx()]
            best_path = [curr.GetIdx()]
            # Look at neighbors that are carbons and allowed and that are not the parent.
            for nbr in curr.GetNeighbors():
                if nbr.GetIdx() == parent:
                    continue
                if nbr.GetAtomicNum() == 6 and allowed_carbon(nbr):
                    path = get_longest_chain(nbr, curr.GetIdx())
                    if len(path) + 1 > len(best_path):
                        best_path = [curr.GetIdx()] + path
            return best_path
        
        chain_path = get_longest_chain(chain_start, acid_c.GetIdx())
        # Main chain is the acid carbon plus the following chain.
        main_chain_idxs = [acid_c.GetIdx()] + chain_path
        chain_length = len(main_chain_idxs)
        
        # Form the fatty acid fragment:
        # Start with the main chain and include the oxygens directly attached to the acid carbon that are part of the acid group.
        fragment_idxs = set(main_chain_idxs)
        for nbr in acid_c.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                fragment_idxs.add(nbr.GetIdx())
        
        frag_size = len(fragment_idxs)
        purity = frag_size / total_heavy if total_heavy > 0 else 0
        
        # Additionally, check for rings: if there is any ring of size >3 that includes atoms not in the fragment then reject.
        extra_ring = False
        for ring in ring_info:
            if len(ring) > 3:
                # If any atom in this ring is not in fragment_idxs, then the molecule has a non‐fatty‐acid ring.
                if not set(ring).issubset(fragment_idxs):
                    extra_ring = True
                    break
        
        # Also measure how many heavy atoms are not in the fragment.
        extra_atoms = total_heavy - frag_size
        
        # Now apply the decision criteria:
        # – The chain length (number of carbons in the main chain including the acid carbon) must exceed 27.
        # – And (either the fatty acid fragment accounts for at least 50% of heavy atoms OR very few extra heavy atoms exist).
        # – And the molecule must not contain any large ring systems outside the fragment.
        if chain_length <= 27:
            reasons.append(f"Chain length is {chain_length} carbons, which does not exceed the C27 threshold")
            continue
        if extra_ring:
            reasons.append("Molecule contains ring(s) larger than 3 atoms outside the fatty acid fragment")
            continue
        if purity < 0.50 and extra_atoms > 5:
            reasons.append(f"Fatty acyl fragment only accounts for {purity:.0%} of heavy atoms (with {extra_atoms} extra atoms)")
            continue
        
        # If all criteria are met:
        return True, f"Chain length is {chain_length} carbons (acid + chain) and fatty acid fragment accounts for {purity:.0%} of heavy atoms."
    
    # If no acid group gives a passing candidate:
    return False, " ; ".join(reasons) if reasons else "No qualifying fatty acid fragment found"

# Example usage:
if __name__ == "__main__":
    # Test with a known ultra-long chain fatty acid (e.g. dotriacontanoic acid: 32 carbons)
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"
    result, explanation = is_ultra_long_chain_fatty_acid(test_smiles)
    print(result, explanation)
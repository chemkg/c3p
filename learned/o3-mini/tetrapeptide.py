"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino–acid residues connected
            by three successive peptide bonds (i.e. a continuous backbone).
            
This implementation does not use a fixed SMARTS pattern for the peptide bond;
instead, it iterates over all bonds in the molecule and looks for bonds where one end is a carbon 
that is “carbonyl‐like” (i.e. is double–bonded to oxygen) and the other is a nitrogen. For each such candidate,
we attempt to assign an “alpha–carbon” (a carbon neighbor of the amide nitrogen other than the carbonyl carbon).
Then we build connectivity between candidates by requiring that the alpha–carbon of one candidate is directly bonded
to the carbonyl carbon of another candidate. Finally, a DFS is used to see if we can “chain” 3 peptide bonds.
"""

from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if the given molecule (as a SMILES string) is a tetrapeptide.
    Strategy:
      1. For every bond in the molecule, check if one atom is a carbon that is carbonyl-like 
         (i.e. has an oxygen double-bonded to it) and the other atom is a nitrogen.
      2. For a candidate peptide bond (i.e. C(=O)-N), attempt to identify the carbon
         attached to the nitrogen that is not the carbonyl carbon (“alpha–carbon”).
      3. Build a connectivity graph between peptide bond candidates: we draw an edge from candidate i
         to candidate j if the alpha–carbon from candidate i is directly bonded to candidate j's carbonyl carbon.
      4. Use depth–first search (DFS) over the candidates graph to see if there is a chain of 3 peptide bonds,
         which would link 4 amino–acid residues.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule appears to be a tetrapeptide, False otherwise.
      str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    peptide_candidates = []  # Each candidate will be a dict with keys: 'carbonyl', 'n', 'alpha'
    
    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # We consider both orders; we want (carbonyl C) bonded to (amide N)
        for c_atom, n_atom in [(a1, a2), (a2, a1)]:
            # Check that c_atom is carbon and n_atom is nitrogen
            if c_atom.GetAtomicNum() != 6 or n_atom.GetAtomicNum() != 7:
                continue
            # Check if c_atom has a double-bonded oxygen (i.e. is a carbonyl carbon)
            has_carbonyl = False
            for nb in c_atom.GetNeighbors():
                if nb.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nb.GetIdx())
                    if bond_to_o is not None and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        has_carbonyl = True
                        break
            if not has_carbonyl:
                continue  # Not a carbonyl carbon
            
            # We now have a candidate carbonyl->nitrogen bond.
            # Attempt to identify an alpha–carbon on the side of the nitrogen.
            alpha = None
            for nb in n_atom.GetNeighbors():
                if nb.GetIdx() == c_atom.GetIdx():
                    continue
                if nb.GetAtomicNum() == 6:  # A carbon neighbor might be the alpha–carbon
                    alpha = nb.GetIdx()
                    break
            if alpha is None:
                continue  # Could not identify an alpha–carbon; skip candidate
            
            peptide_candidates.append({
                'carbonyl': c_atom.GetIdx(),  # carbonyl carbon index
                'n': n_atom.GetIdx(),         # amide nitrogen index
                'alpha': alpha                # chosen alpha–carbon index (on the C-terminal residue)
            })
    
    if not peptide_candidates:
        return False, "No candidate peptide bonds detected."
    
    # Build an adjacency mapping between candidates.
    # For candidates i and j, we draw an edge i -> j if the alpha–carbon of candidate i
    # is directly bonded to the carbonyl carbon of candidate j.
    adjacency = {i: [] for i in range(len(peptide_candidates))}
    for i, cand_i in enumerate(peptide_candidates):
        for j, cand_j in enumerate(peptide_candidates):
            if i == j:
                continue
            # Check if the alpha–carbon of candidate i is directly bonded to candidate j's carbonyl carbon.
            if mol.GetBondBetweenAtoms(cand_i['alpha'], cand_j['carbonyl']) is not None:
                adjacency[i].append(j)
    
    # Use depth-first search (DFS) to look for a chain of 3 peptide bonds.
    # A chain of 3 bonds connects 4 residues.
    def dfs(current, depth, visited):
        # If we have chained 3 bonds, we have found a tetrapeptide backbone.
        if depth == 3:
            return True
        for neighbor in adjacency[current]:
            if neighbor in visited:
                continue
            if dfs(neighbor, depth + 1, visited | {neighbor}):
                return True
        return False
    
    for i in range(len(peptide_candidates)):
        if dfs(i, 1, {i}):
            return True, "Molecule contains a continuous backbone chain with 3 peptide bonds (linking 4 residues) consistent with a tetrapeptide."
    
    return False, "Did not find a continuous backbone chain of 3 peptide bonds connecting 4 residues; peptide bond connectivity does not match a tetrapeptide."

# Example usage (for testing; can be removed in production):
if __name__ == "__main__":
    # Glu-Lys-Trp-Ala example tetrapeptide:
    sample_smiles = "C[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](N)CCC(O)=O)C(O)=O"
    result, reason = is_tetrapeptide(sample_smiles)
    print(result, reason)
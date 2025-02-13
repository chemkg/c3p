"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
Definition: A carbohydrate derivative in which one or more of the oxygens or hydroxy groups 
of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.
Improved heuristic: 
  • For cyclic sugars, search for 5- or 6-membered rings that are non-aromatic and that have 
    (i) at least two exocyclic hydroxyl (–OH) groups (detected by an oxygen with at least one H 
        and bonded only to a ring atom) and 
    (ii) evidence of a thio substitution (either a ring atom that would normally be oxygen replaced 
         by sulfur OR an exocyclic sulfur bound to a ring carbon).
  • For acyclic carbohydrates (open-chain sugars), search for a contiguous chain of 4–6 sp3 carbons 
    that carries at least two hydroxyl groups and one sulfur substituent.
"""

from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string using an improved heuristic.

    Args:
        smiles (str): A SMILES string of the molecule.
    
    Returns:
        (bool, str): True with an explanation if the molecule is classified as a thiosugar;
                     otherwise, False with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens so that we correctly detect –OH groups.
    mol = Chem.AddHs(mol)
    
    # The molecule must contain at least one sulfur atom; if not, it cannot be a thiosugar.
    if not any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
        return False, "No sulfur atom present, so not a thiosugar"
    
    # Helper function: determine if a given neighbor atom is an exocyclic -OH.
    def is_hydroxyl_substituent(neighbor, parent_idx):
        # Check: the neighbor should be oxygen,
        # has at least one hydrogen, and (ideally) is only connected to its parent.
        if neighbor.GetSymbol() != "O":
            return False
        # Count how many heavy (non-H) atoms neighbor this oxygen has (should be one: the ring atom)
        heavy_neighbors = [n for n in neighbor.GetNeighbors() if n.GetSymbol() != "H"]
        if len(heavy_neighbors) != 1:
            return False
        # Ensure that heavy neighbor is indeed the parent atom
        if heavy_neighbors[0].GetIdx() != parent_idx:
            return False
        # Check at least one hydrogen is attached
        if not any(n.GetSymbol() == "H" for n in neighbor.GetNeighbors()):
            return False
        return True

    # Helper function: check if a substituent is an -SR group
    def is_thio_substituent(neighbor, parent_idx):
        if neighbor.GetSymbol() != "S":
            return False
        # For our purposes, if S is attached to the parent and not heavily substituted (or is linked to a small group) we count it.
        # We require that the neighbor S is not itself in a ring with the parent.
        # (We do not examine full substituent details beyond this simple check.)
        if neighbor.IsInRing():
            return False
        return True

    # Check rings (cyclic sugars) first.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    for ring in rings:
        # Only consider rings of 5 or 6 atoms (typical for furanose/pyranose sugars)
        if len(ring) not in (5, 6):
            continue

        # Exclude rings with any aromatic atoms.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # Determine if this ring has a sugar‐like pattern:
        # Count exocyclic -OH groups on any ring atom.
        oh_count = 0
        exo_S_count = 0

        # Also track if one of the ring heteroatoms is S (indicating substitution replacing an oxygen)
        ring_has_S = False

        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() in (8, 16):
                if atom.GetSymbol() == "S":
                    ring_has_S = True
            # Look at neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if is_hydroxyl_substituent(nbr, idx):
                    oh_count += 1
                if is_thio_substituent(nbr, idx):
                    exo_S_count += 1

        # To have a sugar-like profile in these rings we require a reasonable number of -OH groups.
        # Most sugar rings have at least 2 (often 3 or more) exocyclic -OH’s.
        if oh_count < 2:
            continue

        # Now, to be considered a thiosugar, either a ring atom must be sulfur (replacing the expected ring oxygen)
        # or at least one carbon in the ring must carry an S substituent.
        if ring_has_S or exo_S_count >= 1:
            return True, "Thiosugar identified: sugar ring with evidence of a sulfur substitution (in-ring or exocyclic) and sufficient hydroxylation."
    
    # If no cyclic candidate is found, look for acyclic sugar-like patterns.
    # We search for contiguous non‐ring chains of sp3 carbons (length typical for sugars: 4–6) that have several –OH’s.
    atoms = mol.GetAtoms()
    visited = set()
    for atom in atoms:
        # Consider only sp3 carbon atoms that are not in any ring.
        if atom.GetAtomicNum() != 6 or atom.IsInRing():
            continue
        if atom.GetHybridization().name != "SP3":
            continue
        idx = atom.GetIdx()
        if idx in visited:
            continue
        
        # Try to build a chain by a simple DFS (only following sp3 carbons not in rings).
        chain = []
        def dfs(curr, chain_set):
            chain_set.add(curr.GetIdx())
            chain.append(curr)
            for nbr in curr.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()) and nbr.GetHybridization().name=="SP3":
                    if nbr.GetIdx() not in chain_set:
                        dfs(nbr, chain_set)
        dfs(atom, set())
        
        # Only consider chain lengths typical for acyclic sugar fragments.
        if not (4 <= len(chain) <= 6):
            continue
        visited.update(a.GetIdx() for a in chain)

        chain_oh = 0
        chain_S = 0
        for ca in chain:
            for nbr in ca.GetNeighbors():
                if nbr.GetIdx() in [a.GetIdx() for a in chain]:
                    continue
                if is_hydroxyl_substituent(nbr, ca.GetIdx()):
                    chain_oh += 1
                if is_thio_substituent(nbr, ca.GetIdx()):
                    chain_S += 1

        if chain_oh >= 2 and chain_S >= 1:
            return True, "Thiosugar identified: acyclic carbohydrate derivative with sufficient hydroxylation and a sulfur substituent."
    
    return False, "No thiosugar substructure found"

# Example usage:
if __name__ == '__main__':
    # Test with one known thiosugar example, e.g. desulfosinigrin.
    test_smiles = "S([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)C(=NO)CC=C"
    result, reason = is_thiosugar(test_smiles)
    print(result, "->", reason)
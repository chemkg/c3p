"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugate

Definition:
  Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge.
  A bile acid (cholanoic acid) is defined by its tetracyclic (4–ring) fused steroid nucleus (typically three six‐membered rings and one five‐membered ring)
  plus a side chain that has been conjugated (e.g. glycine, taurine, sulfate, glucuronate, sugars, coenzyme A, etc).

This implementation uses an improved heuristic:
  1. Identify a contiguous (fused) set of 4 rings. We require that exactly 4 rings form a single connected component
     (with at least two atoms in common for rings to be “fused”). Also, we check that among these rings there is at least one 5–membered ring and three 6–membered rings.
  2. Look for known conjugation groups (using SMARTS) that attach directly (via at least one bond) to an atom in the fused steroid nucleus.
  
If both features are found, we return True with a reason; otherwise, False plus an explanation.
"""

from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a bile acid conjugate.
    It uses an improved heuristic to:
      1. Identify a fused steroid nucleus (a contiguous set of exactly 4 rings; expected: 1 five‐membered and 3 six‐membered)
      2. Confirm that at least one conjugation group is attached to that nucleus.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an acceptable bile acid conjugate, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- STEP 1: Identify a fused steroid nucleus in the molecule ---
    # Get all rings from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = list(ring_info.AtomRings())
    if not all_rings:
        return False, "No rings detected"
    
    # Build a “ring connectivity graph” where each node is a ring (given as a tuple of atom indices)
    # and an edge exists between two rings if they share 2 or more atoms (i.e. are fused).
    n = len(all_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # If rings i and j share at least 2 atoms, consider them fused.
            if len(set(all_rings[i]).intersection(all_rings[j])) >= 2:
                graph[i].add(j)
                graph[j].add(i)
    
    # Find connected components (clusters of rings)
    def dfs(start, visited):
        stack = [start]
        comp = set()
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            for nb in graph[node]:
                if nb not in visited:
                    stack.append(nb)
        return comp

    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            comp = dfs(i, visited)
            components.append(comp)
    
    # We require that one connected component has exactly 4 fused rings.
    best_comp = None
    for comp in components:
        if len(comp) == 4:
            best_comp = comp
            break
    if best_comp is None:
        return False, "No fused tetracyclic (4-ring) steroid nucleus detected"
    
    # Get the set of atoms in the fused nucleus
    fused_atoms = set()
    for idx in best_comp:
        fused_atoms.update(all_rings[idx])
    
    # Check ring sizes within the fused component
    count_5 = 0
    count_6 = 0
    for idx in best_comp:
        ring_size = len(all_rings[idx])
        if ring_size == 5:
            count_5 += 1
        elif ring_size == 6:
            count_6 += 1
    if count_5 < 1 or count_6 < 3:
        return False, ("Fused nucleus detected but ring sizes do not match expected steroid pattern "
                       "(need at least 1 five-membered and 3 six-membered rings; got: {} 5-membered and {} 6-membered)".format(count_5, count_6))
    
    # --- STEP 2: Check for conjugation groups that are attached to the fused nucleus ---
    # Define SMARTS patterns for groups used in bile acid conjugation.
    conjugation_smarts = {
        "generic amide (e.g. amino acid conjugation)": "C(=O)N",
        "taurine conjugation (anionic)": "NCCS(=O)(=O)[O-]",
        "taurine conjugation (neutral)": "NCCS(=O)(=O)O",
        "sulfate conjugation": "S(=O)(=O)[O-]",
        "glucuronate": "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-]1",
        "sugar (e.g. glucose)": "OC1OC(O)C(O)C(O)C1O"
        # More patterns (e.g. for coenzyme A) can be added as needed.
    }
    
    attached_patterns = []
    for desc, smart in conjugation_smarts.items():
        patt = Chem.MolFromSmarts(smart)
        if patt is None:
            continue
        # Get all matches of this pattern.
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            # For each matched substructure, check whether at least one atom is adjacent to
            # an atom that is part of the fused steroid nucleus.
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                for nbr in atom.GetNeighbors():
                    # If neighbor is in the fused nucleus, we assume attachment.
                    if nbr.GetIdx() in fused_atoms:
                        attached_patterns.append(desc)
                        # Found one attached instance for this pattern, no need to check further.
                        break
                else:
                    continue
                break  # Break out of the match loop if found
    if not attached_patterns:
        return False, "Fused steroid nucleus detected but no recognized conjugation group is attached to it"
    
    reason = ("Fused tetracyclic steroid nucleus detected (1 five-membered and 3 six-membered rings) "
              "with conjugation group(s): " + ", ".join(sorted(set(attached_patterns))))
    return True, reason


# For testing purposes, here are some example SMILES strings.
if __name__ == "__main__":
    smiles_examples = [
        # taurocholic acid (should be True)
        "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCC(=O)NCCS(O)(=O)=O",
        # tauro-omega-muricholic acid (should be True)
        "[C@]12([C@]([C@]3([C@]([C@@H]([C@@H](CCC(NCCS(=O)(=O)O)=O)C)CC3)(C)CC1)[H])([C@H]([C@H](O)[C@]4([C@@]2(CC[C@H](C4)O)C)[H])O)[H])[H]",
        # a false positive case from previous attempt (Epicoccin U-like structure):
        "S1S[C@]23N([C@@H]4[C@@H](O)[C@@H](O)CC([C@@H]4C2)=O)C([C@]15N([C@@H]6[C@@H](OC)[C@@H](O)CC([C@@H]6C5)=O)C3=O)=O"
    ]
    for s in smiles_examples:
        result, reason = is_bile_acid_conjugate(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("----------")
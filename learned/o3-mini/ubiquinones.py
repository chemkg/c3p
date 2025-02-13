"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: Ubiquinones – any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone.
Ubiquinones usually carry a polyprenyl side chain (an isoprenoid chain) at position 6.
This improved version:
  • Iterates over each six-membered ring to detect a quinone core by:
      – Detecting exactly two ring carbons with exocyclic carbonyls.
      – Ensuring that among the remaining ring atoms there are at least two oxygen substituents
        (indicative of a dimethoxy/hydroxy pattern) and at least one methyl substituent.
  • Then it attempts to identify a unique external branch (the candidate position-6 side chain)
    that is not already a carbonyl/oxygen/methyl substituent.
  • It “grows” the branch (using a DFS while avoiding the core) and requires that some double bond is present
    (to reject saturated chains). It further counts prenyl (isoprene) repeating units using a SMARTS match.
  • If the candidate core and prenyl branch are found, the function returns True along with a message
    stating the number of isoprenoid unit(s) detected.
If no candidate ring passes these tests the function returns False with an explanation.
"""

from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule belongs to the class of ubiquinones based on its SMILES string.
    It looks for a 2,3-dimethoxy-5-methylbenzoquinone-derived core (a six-membered ring with two carbonyl groups,
    at least two oxy functions on non-carbonyl atoms, and at least one methyl substituent) and requires that a polyprenyl 
    (isoprenoid) branch attached (at position 6) is present. The branch must contain at least one double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ubiquinone, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS query for one prenyl (isoprene) unit: a fragment with an alkene and a methyl branch.
    prenyl_query = Chem.MolFromSmarts("C/C=C(/C)")
    
    # Helper to check for a double-bonded oxygen (i.e. a carbonyl bond)
    def is_carbonyl(bond, nbr):
        return (nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    
    # Helper to decide if an oxygen substituent is methoxy or hydroxy.
    def is_oxy_substituent(oatom, from_atom_idx):
        # Expects an oxygen not connected via a double bond.
        if oatom.GetAtomicNum() != 8:
            return False
        # Either a hydroxy (at least one hydrogen) or methoxy (attached to a CH3)
        if oatom.GetTotalNumHs() >= 1:
            return True
        # Alternatively, check for methoxy: one neighbor besides the core should be a methyl.
        for nb in oatom.GetNeighbors():
            if nb.GetIdx() == from_atom_idx:
                continue
            if nb.GetAtomicNum() == 6 and nb.GetDegree() == 1 and nb.GetTotalNumHs() == 3:
                return True
        return False
    
    # Helper to check if a carbon atom is a methyl group.
    def is_methyl(catom):
        return (catom.GetAtomicNum() == 6 and catom.GetDegree() == 1 and catom.GetTotalNumHs() == 3)
    
    # For a given branch attachment, collect all atom indices reachable off the core (avoid core atoms).
    def get_branch_atom_indices(start_idx, core_set):
        branch_atoms = set()
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in branch_atoms:
                continue
            branch_atoms.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_set:
                    continue
                if nbr.GetIdx() not in branch_atoms:
                    stack.append(nbr.GetIdx())
        return branch_atoms
    
    # Check if the branch (given by a set of atom indices) contains at least one double bond.
    def branch_has_double_bond(branch_atom_indices):
        for idx in branch_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in branch_atom_indices:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    return True
        return False

    # Our strategy:
    # Iterate over each 6-membered ring in the molecule. For each, we try to detect:
    #   1. Exactly two ring carbons with an exocyclic double-bonded oxygen (carbonyl groups)
    #   2. Among remaining ring carbons, at least 2 oxygen substituents (methoxy/hydroxy) and at least 1 methyl substituent.
    #   3. A unique external branch (other than those already counted) which qualifies as a polyprenyl side chain (with unsaturation).
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        if len(ring) != 6:
            continue  # We only consider six-membered rings.
        ring_set = set(ring)
        
        carbonyl_indices = set()  # ring atom indices that have an exocyclic C=O.
        oxy_substituent_count = 0
        methyl_substituent_count = 0
        # We will record, for each ring atom, the external neighbor indices and classify them.
        branch_candidates = set()  # potential attachment atoms for the prenyl chain.
        
        # First, mark carbonyls in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only carbon atoms considered for the core.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if is_carbonyl(bond, nbr):
                    carbonyl_indices.add(idx)
                    break  # Only count once per ring atom.
        
        # We require exactly two ring carbons with carbonyl groups.
        if len(carbonyl_indices) != 2:
            continue
        
        # Now, on each ring atom, examine external substituents.
        # For atoms not already classified as carbonyl-bearing, record oxygen and methyl substituents.
        # Also record any external carbon that is not a simple methyl (could be prenyl branch).
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Skip if this ring atom already gave a carbonyl bond.
            # (Assume that if an atom is carbonyl-bearing, its substituent is "taken".)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                # Check if this neighbor was already accounted for by being a carbonyl (by checking bond type)
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    continue  # already in carbonyl
                # Classify oxygens:
                if nbr.GetAtomicNum() == 8:
                    if is_oxy_substituent(nbr, atom.GetIdx()):
                        oxy_substituent_count += 1
                    continue
                # Classify methyl substituents:
                if nbr.GetAtomicNum() == 6:
                    if is_methyl(nbr):
                        methyl_substituent_count += 1
                    else:
                        # Candidate branch for prenyl side chain.
                        branch_candidates.add(nbr.GetIdx())
        
        # For the core matching we require at least 2 oxygen substituents and at least 1 methyl group.
        if oxy_substituent_count < 2 or methyl_substituent_count < 1:
            continue
        
        # Next, from branch candidates, we try to pick exactly one candidate branch that is adjacent to the core.
        # (If several exist, we will process each in turn.)
        for branch_start in branch_candidates:
            branch_atoms = get_branch_atom_indices(branch_start, ring_set)
            # Check that the branch has at least one double bond (to distinguish unsaturated polyprenyl from a saturated alkyl chain).
            if not branch_has_double_bond(branch_atoms):
                continue
            # Now, attempt to count the isoprene (prenyl) units in the branch.
            # We extract a sub-molecule corresponding to the branch.
            branch_atom_list = list(branch_atoms)
            try:
                branch_mol = Chem.PathToSubmol(mol, branch_atom_list)
            except Exception:
                continue
            prenyl_matches = branch_mol.GetSubstructMatches(prenyl_query)
            n_prenyl = len(prenyl_matches)
            if n_prenyl < 1:
                continue  # No prenyl unit detected on this branch.
            # If we get here, we have a candidate quinone core with a prenyl branch.
            reason = (f"Contains 2,3-dimethoxy-5-methylbenzoquinone-derived core with "
                      f"{n_prenyl} isoprenoid unit(s) detected")
            return True, reason

    # If no candidate ring passed our screening, then reject.
    return False, ("Does not contain a suitable 2,3-dimethoxy-5-methylbenzoquinone derived core "
                   "with proper substituents and a prenyl side chain")

# Example usage (uncomment the block below to test):
# if __name__ == "__main__":
#     tests = [
#         # True positives:
#         "O=C1C(OC)=C(OC)C(=O)C(C\\C=C(/CC/C=C(/C[C@H](O)C=C(C)C)\\C)\\C)C1=O",  # Antroquinonol N
#         "COC1=C(OC)C(=O)C(CC=C(C)C)=C(C)C1=O",                                # ubiquinone-1
#         # False negatives (expected ubiquinones but our improved method might need further tuning):
#         "S(C=1C(=O)C(OC)=C(OC)C(C1CC(=O)OC)=O)C",                              # Coptirhoquinone A
#         "O(C=1C(=O)C([C@@H](C2=CC=C(O)C=C2)C=C)=CC(=O)C1OC)C",                  # (R)-4'-Hydroxy-3,4-dimethoxydalbergione
#     ]
#     for s in tests:
#         result, explanation = is_ubiquinones(s)
#         print(result, explanation)
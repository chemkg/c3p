"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: Catecholamine (4-(2-Aminoethyl)pyrocatechol and substitution derivatives)

A catecholamine is defined here as a molecule that contains an aromatic benzene ring (a 6-membered aromatic ring)
that carries a catechol-like motif – that is, two adjacent substituents that are oxygen atoms (they can be –OH or 
derivatized such as –OSO3H, –OCH3, etc.), plus an attached ethylamine chain (or a close analogue) originating from one of 
the ring atoms.
"""

from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine derivative 
    (e.g., dopamine, adrenaline, dobutamine, etc.) based on its SMILES string.
    It does so by detecting two local features:
      1. The presence of a benzene (aromatic 6-membered) ring in which two adjacent ring atoms have a substituent
         that is an oxygen atom (the catechol motif). The oxygen substituents can be modified (e.g. sulfate, methoxy, etc.).
      2. The presence of an attached ethylamine-like chain; that is a short chain of 2–3 bonds from a ring atom leading to a nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        tuple(bool, str): True with a classification reason if the molecule qualifies as a catecholamine derivative;
                          False with a reason otherwise.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule sanitization failed: " + str(e)
    
    # Get ring information. We are interested in aromatic 6-membered rings.
    ring_info = mol.GetRingInfo().AtomRings()
    
    catechol_found = False
    chain_found = False
    
    # Helper function: from a starting atom (which is an exocyclic neighbor to a ring atom), 
    # search for a path (of length at most max_depth) that leads to a nitrogen atom.
    # We require that at least 2 bonds (non-zero depth) are traversed.
    def has_ethylamine_chain(start_atom, max_depth=3):
        from collections import deque
        # each item: (atom, depth, came_from_idx)
        queue = deque()
        queue.append((start_atom, 1, None))
        while queue:
            atom, depth, prev = queue.popleft()
            # if we have traversed at least two bonds and find a nitrogen, we say found chain.
            if depth >= 2 and atom.GetAtomicNum() == 7:
                return True
            if depth < max_depth:
                for nbr in atom.GetNeighbors():
                    # avoid going back along the bond we came from
                    if nbr.GetIdx() == prev:
                        continue
                    # We restrict the search to aliphatic atoms (or at least not aromatic since chain atoms are usually not aromatic)
                    # However, we allow heteroatoms in the chain if needed.
                    # We do not traverse into rings (to avoid accidental ring-closing to an aromatic ring).
                    if nbr.GetIdx() in current_ring:
                        continue
                    queue.append((nbr, depth+1, atom.GetIdx()))
        return False

    # Loop over aromatic rings of size 6.
    for ring in ring_info:
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # First, check for catechol motif: at least one pair of adjacent ring atoms (in the ring order) 
        # must have an exocyclic oxygen substituent.
        adjacent_oxy = False
        # record for later chain search: the ring atoms that have an oxygen substituent.
        oxy_ring_atoms = set()
        # iterate over atoms in the ring
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # check each neighbor that is not in the ring:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # if neighbor is oxygen (atomic num 8), we count it.
                if nbr.GetAtomicNum() == 8:
                    oxy_ring_atoms.add(idx)
                    break  # one oxygen substituent is enough
        
        # Now check if there exists a pair of adjacent ring atoms (cyclically adjacent in the ring) 
        # that both have an oxygen substituent.
        ring_list = list(ring)
        n_atoms = len(ring_list)
        for i in range(n_atoms):
            a1 = ring_list[i]
            a2 = ring_list[(i+1) % n_atoms]  # adjacent in the cycle
            if a1 in oxy_ring_atoms and a2 in oxy_ring_atoms:
                adjacent_oxy = True
                break
        
        if not adjacent_oxy:
            continue  # this ring does not have the catechol motif
        
        # Mark that we have a catechol (or catechol-derivative) ring.
        catechol_found = True
        
        # Save this ring as 'current_ring' for chain search (so that we do not traverse into other ring atoms)
        current_ring = set(ring)
        
        # Next, look for an ethylamine chain attached to this ring.
        # For each ring atom, examine its exocyclic neighbors.
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in current_ring:
                    continue  # ignore neighbors inside the ring
                # we require that the chain starts from a carbon (usually)
                if nbr.GetAtomicNum() != 6:
                    continue
                # We also prefer that this first neighbor is not aromatic (i.e. an aliphatic carbon).
                if nbr.GetIsAromatic():
                    continue
                # Now, perform a limited-depth search from this neighbor to see if we can hit a nitrogen (chain end).
                if has_ethylamine_chain(nbr, max_depth=3):
                    chain_found = True
                    break
            if chain_found:
                break
        
        # If for this ring we found both the catechol motif and an attached ethylamine chain, classify as catecholamine.
        if catechol_found and chain_found:
            return True, "Molecule contains a catechol aromatic ring (adjacent oxygen substituents) with an attached short chain leading to a nitrogen"
    
    # If no ring satisfied both conditions, construct a reason
    reasons = []
    if not catechol_found:
        reasons.append("No aromatic six‐membered ring with two adjacent oxygen‐bearing substituents (catechol motif) was found")
    elif not chain_found:
        reasons.append("Catechol ring found but no attached ethylamine (chain that reaches a nitrogen) was found")
    return False, "; ".join(reasons)
    
# Example usage (for testing):
if __name__ == '__main__':
    # List some test SMILES strings with names for debugging
    test_molecules = [
        ("(S)-dobutamine", "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1"),
        ("Epinephrine sulfate", "S(OC1=C(O)C=C([C@@H](O)CNC)C=C1)(O)(=O)=O"),
        ("Normetanephrine", "C=1(C=C(C(=CC1)O)OC)C(O)CN"),
        ("Dopamine", "NCCc1ccc(O)c(O)c1"),
    ]
    for name, smi in test_molecules:
        result, reason = is_catecholamine(smi)
        print(f"{name}: {result} - {reason}")
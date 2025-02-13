"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: Catecholamine (4-(2-Aminoethyl)pyrocatechol and substitution derivatives)

A catecholamine is defined here as a molecule containing a benzene (aromatic 6-membered) ring with two adjacent hydroxyl (–OH) groups 
(i.e. a catechol moiety) AND an exocyclic ethylamine chain (–CH2–CH2–NHx) attached to one of the ring atoms.
"""

from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine derivative (e.g. dopamine, adrenaline, dobutamine, etc.) 
    based on its SMILES string. It does so by checking for the presence of a benzene ring that has 
    adjacent hydroxyl substituents (a catechol group) and an ethylamine chain attached to that ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True and a reason if the molecule qualifies as a catecholamine derivative;
                     False and a reason otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule (adjust hydrogens, valences, etc.)
    Chem.SanitizeMol(mol)
    
    ring_info = mol.GetRingInfo().AtomRings()
    catechol_found = False
    chain_found = False

    # For each ring in the molecule
    for ring in ring_info:
        # We only consider rings of size 6
        if len(ring) != 6:
            continue

        # Check that the ring is aromatic
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # For each atom in the ring, mark if it has an -OH substituent
        # (we decide an -OH is an oxygen attached to the ring atom that has at least one implicit hydrogen)
        oh_atoms = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider substituents not in the ring
                if nbr.GetIdx() in ring:
                    continue
                # Check if this neighbor is oxygen and has at least one hydrogen (using total H count)
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    oh_atoms.add(idx)
                    break  # one -OH is enough

        # Now check if at least two adjacent ring atoms (neighbors in the ring cycle) both have -OH.
        adjacent_oh = False
        ring_list = list(ring)  # ordered list of indices for the ring atoms
        n_atoms = len(ring_list)
        for i in range(n_atoms):
            curr_idx = ring_list[i]
            next_idx = ring_list[(i + 1) % n_atoms]  # cyclic adjacent
            if curr_idx in oh_atoms and next_idx in oh_atoms:
                adjacent_oh = True
                break

        if not adjacent_oh:
            continue  # this ring doesn't appear to have a catechol motif
        
        # Mark that we have found a catechol ring.
        catechol_found = True

        # Now check for the side chain: an ethylamine chain attached to any ring atom.
        # We search for a branch: RingAtom -- CH2 (aliphatic carbon) -- CH2 (aliphatic carbon) -- N (any amino)
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # only consider exocyclic substituents
                # The first chain atom must be a carbon and not aromatic.
                if nbr.GetAtomicNum() != 6 or nbr.GetIsAromatic():
                    continue
                chain1 = nbr
                # Look for a second chain carbon attached to chain1 (not going back to ring_atom)
                for nbr2 in chain1.GetNeighbors():
                    if nbr2.GetIdx() == ring_atom.GetIdx():
                        continue  # do not backtrack
                    if nbr2.GetAtomicNum() != 6 or nbr2.GetIsAromatic():
                        continue
                    chain2 = nbr2
                    # Now, from chain2, check for a neighbor that is a nitrogen (any substitution is allowed).
                    for nbr3 in chain2.GetNeighbors():
                        # Exclude backtracking to chain1.
                        if nbr3.GetIdx() == chain1.GetIdx():
                            continue
                        if nbr3.GetAtomicNum() == 7:
                            chain_found = True
                            break
                    if chain_found:
                        break
                if chain_found:
                    break
            if chain_found:
                break

        # If we found both motifs in this ring, we classify as catecholamine.
        if catechol_found and chain_found:
            return True, "Molecule contains a catechol ring (adjacent hydroxyl groups) with an attached ethylamine chain"

    # If none of the rings satisfied both conditions:
    reasons = []
    if not catechol_found:
        reasons.append("No aromatic six-membered ring with adjacent hydroxyl groups (catechol group) was found")
    elif not chain_found:
        reasons.append("No attached ethylamine (–CH2–CH2–NHx) chain was found on the catechol ring")
        
    reason_str = "; ".join(reasons)
    return False, reason_str
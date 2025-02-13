"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: Hydroxynaphthoquinone
Definition:
  A hydroxynaphthoquinone is any naphthoquinone in which the naphthoquinone moiety 
  (the bicyclic core originally derived from naphthalene) has at least one hydroxy (-OH) group attached.
  A naphthoquinone is defined as a fused bicyclic ring system (two six‐membered rings sharing two atoms)
  that has at least two carbonyl (C=O) groups attached (where the oxygen is not part of the core).

Approach:
  1. Parse the SMILES string and add explicit hydrogens (so –OH groups become visible).
  2. Create a “core” molecule with hydrogens removed (to allow proper ring detection).
  3. Examine all rings in the core to find two six‐membered rings sharing exactly two atoms.
     These combined atoms are taken as the candidate naphthalene core.
  4. For each carbon in the candidate core (from the original mol with H’s),
     check for attached carbonyl groups (a double bond from the core carbon to an oxygen not in the core).
  5. Also check for at least one hydroxy substituent: an oxygen (bonded by a single bond) that has at least one hydrogen.
  6. Return True if there are at least 2 carbonyls and at least one hydroxy substituent on the core.
  
  If any of these tests fail (or the core cannot be found), return False with an explanation.
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines whether a given molecule qualifies as a hydroxynaphthoquinone.

    A hydroxynaphthoquinone is defined as a naphthoquinone (a fused bicyclic six-membered ring system
    with at least two carbonyl (C=O) groups attached – with the oxygens outside the ring core) that is substituted
    by at least one hydroxy (-OH) group on the core.

    Args:
        smiles (str): The SMILES string for the molecule.

    Returns:
        bool: True if the molecule qualifies as a hydroxynaphthoquinone, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to properly see -OH substituents.
    mol = Chem.AddHs(mol)
    
    # For ring analysis, remove explicit hydrogens so that aromaticity and ring info is standard.
    core_mol = Chem.RemoveHs(mol)
    ring_info = core_mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Find candidate naphthalene core: two 6-membered rings sharing exactly 2 atoms.
    candidate_core = None
    six_membered = [set(ring) for ring in atom_rings if len(ring) == 6]
    for i in range(len(six_membered)):
        for j in range(i+1, len(six_membered)):
            # Check if these two rings share exactly 2 atoms.
            intersection = six_membered[i].intersection(six_membered[j])
            if len(intersection) == 2:
                # Candidate core is the union of the atoms in the two rings.
                candidate_core = six_membered[i].union(six_membered[j])
                break
        if candidate_core is not None:
            break

    if candidate_core is None:
        return False, "No fused bicyclic six-membered ring system (naphthalene core) found"

    # Now count carbonyl (C=O) groups attached to atoms in the candidate core.
    carbonyl_count = 0
    hydroxy_count = 0
    # To ensure we count one substitution per core atom, we record indices that already have a substitution.
    atoms_with_oh = set()
    atoms_with_co = set()
    
    for idx in candidate_core:
        atom = mol.GetAtomWithIdx(idx)
        # We check only carbon atoms in the core.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Check each bond from this core atom.
        for bond in atom.GetBonds():
            # Get the neighboring atom.
            nbr = bond.GetOtherAtom(atom)
            nbr_idx = nbr.GetIdx()
            # Consider only bonds from the core to outside the core.
            if nbr_idx in candidate_core:
                continue
            
            # Check for a carbonyl group: a double bond from C to an oxygen.
            if bond.GetBondType() == Chem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                # It is a carbonyl if the oxygen is not part of the core.
                atoms_with_co.add(idx)
            # Check for hydroxy: an oxygen singly bonded that carries at least one hydrogen.
            if bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetAtomicNum() == 8:
                # Sometimes an -OH group might be depicted with explicit hydrogens.
                # Check for at least one attached hydrogen.
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 1:  # hydrogen
                        atoms_with_oh.add(idx)
                        break

    carbonyl_count = len(atoms_with_co)
    hydroxy_count = len(atoms_with_oh)
    
    # According to our definition of naphthoquinone, we require at least 2 carbonyl groups.
    if carbonyl_count < 2:
        return False, ("Naphthalene core found but has only {} carbonyl group(s) attached; "
                       "need at least 2 to be a naphthoquinone".format(carbonyl_count))
    
    # For a hydroxynaphthoquinone, we require at least one hydroxy group on the core.
    if hydroxy_count < 1:
        return False, ("Naphthoquinone core has {} carbonyl group(s) but lacks a hydroxy (-OH) substituent "
                       "on the core".format(carbonyl_count))
    
    return True, ("Found fused bicyclic core with {} carbonyl group(s) and {} hydroxy substituent(s) on the core, "
                  "qualifying as a hydroxynaphthoquinone".format(carbonyl_count, hydroxy_count))

# Example usage:
# smiles_example = "Oc1ccc(O)c2C(=O)C=CC(=O)c12"  # naphthazarin, for instance
# print(is_hydroxynaphthoquinone(smiles_example))
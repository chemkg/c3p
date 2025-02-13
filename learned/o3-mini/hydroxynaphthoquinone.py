"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: Hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety is substituted by at least one hydroxy group.
The approach is to:
  1. Parse the SMILES and add explicit hydrogens.
  2. Use the molecule’s ring system (via GetRingInfo) to find aromatic 6-membered rings.
  3. Identify candidate naphthalene cores by finding pairs of aromatic 6-membered rings that share at least 2 atoms.
  4. For each candidate core, count the number of carbonyl groups attached (i.e. an aromatic carbon that is double‐bonded
     to an oxygen that is not itself part of the fused ring) and count hydroxy substituents on atoms of the fused system.
  5. Qualify it as a hydroxynaphthoquinone if the candidate core shows at least two carbonyl groups and at least one –OH substituent.
"""

from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines whether a molecule is a hydroxynaphthoquinone based on its SMILES string.

    A hydroxynaphthoquinone is defined as a naphthoquinone (a fused bicyclic aromatic system that has at least
    two carbonyl groups incorporated into the ring) in which at least one carbon of the fused system is substituted by a hydroxy (-OH) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule qualifies as hydroxynaphthoquinone, otherwise False.
        str: Detailed explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that -OH groups become explicit.
    mol = Chem.AddHs(mol)

    ring_info = mol.GetRingInfo()
    if not ring_info:
        return False, "No ring information found in the molecule"

    atom_rings = ring_info.AtomRings()
    # Filter rings: we want 6-membered aromatic rings.
    aromatic_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            # check if all atoms in the ring are aromatic (and are carbons)
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_rings.append(set(ring))
    
    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found; naphthalene core missing"

    candidate_cores = []  # each candidate is a set of atom indices representing a fused pair
    # Look for pairs of aromatic rings that are fused (share at least 2 atoms).
    for i in range(len(aromatic_rings)):
        for j in range(i+1, len(aromatic_rings)):
            intersec = aromatic_rings[i] & aromatic_rings[j]
            if len(intersec) >= 2:
                # The union of the two rings forms a candidate naphthalene core.
                candidate_cores.append(aromatic_rings[i] | aromatic_rings[j])
    
    if not candidate_cores:
        return False, "No fused aromatic ring system (naphthalene core) found"

    # For each candidate naphthalene core, assess for quinone and hydroxy substituents.
    for core in candidate_cores:
        carbonyl_count = 0
        hydroxy_count = 0
        # For each atom in the core
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider aromatic carbons from the core.
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                continue

            # Check bonds of the atom. We consider a bond from an aromatic carbon to an oxygen (outside the core) that is a double bond as a carbonyl.
            for bond in atom.GetBonds():
                # Look for double bonds to oxygen.
                if bond.GetBondTypeAsDouble() == 2.0:
                    nbr = bond.GetOtherAtom(atom)
                    # Ensure the neighbor is O and is not part of the fused core.
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in core:
                        carbonyl_count += 1
                        break  # count each core carbon only once
            # Check for hydroxy substituent. Look at neighbors for an -OH group.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in core:
                    # Check if this oxygen is bonded to at least one hydrogen.
                    # (Since we added H, the explicit hydrogen count should be available)
                    if nbr.GetTotalNumHs() > 0:
                        hydroxy_count += 1
                        break  # count each core carbon only once

        # Now decide based on counts.
        if carbonyl_count < 2:
            # Candidate core does not show a naphthoquinone pattern.
            continue  # try the next candidate core
        if hydroxy_count < 1:
            return False, ("Fused aromatic core with sufficient carbonyl groups detected ({} carbonyl(s)) "
                           "but no hydroxy substituent on the core".format(carbonyl_count))
        # If criteria are met:
        return True, ("Found fused aromatic core (naphthalene-like) with {} carbonyl group(s) and "
                      "{} hydroxy substituent(s) on the core".format(carbonyl_count, hydroxy_count))

    # If none of the candidate cores met the combined criteria:
    return False, "Fused aromatic core found but it does not have both the required quinone (>=2 carbonyls) and hydroxy substituent(s)"

# (The module can be tested by calling the function with one or more SMILES strings)
# Example:
# print(is_hydroxynaphthoquinone("OC1=CC(=O)c2ccccc2C1=O"))
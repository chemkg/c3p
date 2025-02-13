"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: Catechols (compounds containing an o-diphenol component)

A catechol is defined here as any compound that contains a benzene ring
with two hydroxyl groups on two adjacent positions. Prior attempts that relied
on a fixed SMARTS pattern ("[c]([OX2H])c([OX2H])") turned out to be too liberal,
leading to many false positives (and missing some deprotonated or sterically subtle cases).
Here we try a different approach:
  1. Parse the molecule.
  2. Retrieve every aromatic ring that has exactly 6 atoms.
  3. For each such ring, examine each pair of adjacent ring atoms.
     For a given atom we check whether it bears a hydroxyl group – defined here
     as a bonded oxygen that (after implicit H are added) has exactly one neighbor.
  4. If at least one ring has a pair of adjacent -OH groups on a benzene ring,
     we classify it as catechol.

If any severe error is encountered or no good decision can be made,
the function may return (None, None).
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains an o-diphenol (catechol) component based on its SMILES.
    The procedure looks for a six‐membered aromatic ring (benzene) that has two adjacent
    atoms each of which bears a hydroxyl substituent.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a candidate catechol moiety is identified, False if not.
        str: Explanation of the decision.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # (We could add hydrogens explicitly to help with hydroxyl detection.)
    mol = Chem.AddHs(mol)

    # Helper function: given an atom, does it have a hydroxyl (-OH) substituent?
    # We define a hydroxyl group as an oxygen atom (atomic number 8) attached
    # to the candidate atom, where that oxygen is bonded only to that atom (degree == 1)
    # and carries at least one hydrogen.
    def has_hydroxyl(atom):
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check: oxygen attached only to this heavy atom.
                # (In RDKit the number of (implicit+explicit) H's is set on the oxygen.)
                if nbr.GetDegree() == 1 and nbr.GetTotalNumHs() >= 1:
                    return True
        return False

    # Look at every ring in the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    # Loop over each ring
    for ring in atom_rings:
        # Only consider 6-membered rings.
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is aromatic (i.e. a benzene ring)
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # For the atoms in the ring, check each cyclic adjacent pair.
        # Because the ring is cyclic, we loop with wrap‐around.
        n_atoms = len(ring)
        for i in range(n_atoms):
            idx1 = ring[i]
            idx2 = ring[(i + 1) % n_atoms]
            atom1 = mol.GetAtomWithIdx(idx1)
            atom2 = mol.GetAtomWithIdx(idx2)
            if has_hydroxyl(atom1) and has_hydroxyl(atom2):
                return True, "Found catechol moiety on benzene ring with adjacent -OH groups"

    # If we get here, no appropriate ring was found.
    return False, "No o-diphenol (catechol) substructure found"

# The following test block (if uncommented) will run example molecules.
if __name__ == '__main__':
    # Some example SMILES strings from the outcomes provided:
    test_examples = [
        # True positive examples (some true catechol compounds)
        "OC1=C(O)C=CC=C1CCCC/C=C\\C/C=C\\CCCCCCCC2=C(O)C(O)=CC=C2",  # Gerronemin F
        "COc1cc(CCc2ccc(O)c(O)c2)cc(O)c1O",  # dendrocandin E
        "O[C@@H](CC\\C=C\\c1ccccc1)CCc1ccc(O)c(O)c1",  # (-)-(3S)-1-(3,4-dihydroxyphenyl)-7-phenyl-(6E)-6-hepten-3-ol
        # A false positive example from previous attempt (e.g. D-dopa)
        "N[C@H](Cc1ccc(O)c(O)c1)C(O)=O",  # D-dopa
    ]
    for s in test_examples:
        classification, reason = is_catechols(s)
        print("SMILES: {}\nClassification: {}\nReason: {}\n".format(s, classification, reason))
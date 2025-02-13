"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoid – Members of the class of tetrahydrochromenochromene.
Definition: Rotenoids consist of a cis‐fused tetrahydrochromenochromene skeleton
and its substituted derivatives. Our refined heuristic uses ring analysis:
  1. Retrieve all rings in the molecule.
  2. Loop over pairs of fused rings (sharing at least 2 atoms).
  3. For each pair, compute the fraction of atoms in each ring that are aromatic.
  4. Identify the ring with the higher aromatic fraction and require that:
       - Its aromatic fraction is at least 0.75 (chromene‐like).
       - The difference between the aromatic fractions is at least 0.10.
       - The more aromatic ring contains at least one oxygen (atomic number 8)
         that is a member of that ring.
If any ring pair satisfies these conditions, the molecule is classified as a rotenoid.
"""

from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is likely to be a rotenoid based on its SMILES string.
    The heuristic approach is:
      1. Get all rings (lists of atom indices) in the molecule.
      2. For every pair of rings that are fused (i.e., share at least 2 atoms),
         calculate the aromatic fractions.
      3. Identify the more aromatic ring of the pair. It must have an aromatic fraction
         at least 0.75 and differ from its partner by at least 0.10.
      4. In addition, the more aromatic ring must have at least one oxygen atom as a ring member.
    Args:
        smiles (str): SMILES representation of the molecule.
    Returns:
        bool: True if the molecule likely belongs to the rotenoid class, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Retrieve all rings in the molecule as tuples of atom indices.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    aromatic_threshold = 0.75
    aromatic_diff_threshold = 0.10

    # Loop over all pairs of rings to find a fused bicyclic system.
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            # Check if the rings share at least 2 atoms (i.e., fused via a bond).
            shared_atoms = set(rings[i]).intersection(rings[j])
            if len(shared_atoms) < 2:
                continue  # Not fused strongly enough

            # Get atom objects for each ring.
            ring1_atoms = [mol.GetAtomWithIdx(idx) for idx in rings[i]]
            ring2_atoms = [mol.GetAtomWithIdx(idx) for idx in rings[j]]

            # Calculate aromatic fractions for each ring.
            arom1 = sum(1 for atom in ring1_atoms if atom.GetIsAromatic())
            arom2 = sum(1 for atom in ring2_atoms if atom.GetIsAromatic())
            frac1 = arom1 / len(ring1_atoms)
            frac2 = arom2 / len(ring2_atoms)

            # Identify the more aromatic ring and the less aromatic ring.
            if frac1 >= frac2:
                high_frac, low_frac = frac1, frac2
                high_ring_atoms = ring1_atoms
            else:
                high_frac, low_frac = frac2, frac1
                high_ring_atoms = ring2_atoms

            # Check if the more aromatic ring meets the aromaticity threshold 
            # and the difference in aromaticity is sufficiently high.
            if high_frac < aromatic_threshold or (high_frac - low_frac) < aromatic_diff_threshold:
                continue

            # Check that the more aromatic ring contains at least one oxygen atom
            # as an integral member of the ring.
            contains_oxygen = any(atom.GetAtomicNum() == 8 for atom in high_ring_atoms)
            if not contains_oxygen:
                continue

            # If all criteria are met, classify the molecule as a rotenoid.
            return True, ("Molecule contains a fused bicyclic system with one predominantly aromatic ring "
                          "(aromatic fraction {:.2f}) and another less aromatic ring (fraction {:.2f}), "
                          "with at least one oxygen atom in the aromatic ring. This is consistent with a rotenoid skeleton."
                          .format(high_frac, low_frac))
    
    # If no pair of rings passes the criteria, classify as not rotenoid.
    return False, "Molecule does not contain the fused bicyclic system characteristic of rotenoids"

# Example usage (testing a few molecules)
if __name__ == "__main__":
    test_examples = {
        "Millettosin": "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C6=C(OC(C=C6)(C)C)C=C5",
        "LSM-4404": "CC(C)C1CC2=C(O1)C=CC3=C2OC4=C(C3=O)C5=CC(=C(C=C5OC4)OC)OC",
        "Tephrosin": "[H][C@@]12COc3cc(OC)c(OC)cc3[C@]1(O)C(=O)c1ccc3OC(C)(C)C=Cc3c1O2",
        "Stemonone": "O1C2=C(C=3C(OC2=O)=CC(OC)=C(OC)C3)C(=O)C=4C1=CC(OC)=CC4O"
    }
    for name, smiles in test_examples.items():
        result, reason = is_rotenoid(smiles)
        print(f"{name}: {result} -- {reason}")
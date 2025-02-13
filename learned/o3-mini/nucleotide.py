"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This program uses several structural criteria:
1. The molecule must contain at least one phosphorus atom.
2. It must contain a sugar moiety – here we identify a five‐membered furanose ring (4 carbons and 1 oxygen).
3. It must contain an aromatic nucleobase (an aromatic ring with at least 2 nitrogen atoms) that is not the sugar ring.
4. A phosphate must be attached to the sugar. Instead of demanding that the phosphate oxygen be in the furanose ring,
   we now allow the phosphate to be attached to an oxygen that either belongs to the sugar ring or is directly attached to a sugar ring atom.
"""

from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    The algorithm checks for:
      - Presence of phosphorus (P) for the phosphate group.
      - A five-membered sugar ring (furanose) defined as a ring with 5 atoms composed of 4 carbons and 1 oxygen.
      - An aromatic nucleobase ring (not the sugar) containing at least 2 nitrogen atoms.
      - A phosphate group attached (via an oxygen) to the sugar (either directly as part of the sugar ring or
        to an oxygen that is bound to a sugar ring atom, corresponding to condensation at the 3' or 5' OH).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a nucleotide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for phosphorus atoms (indicates potential phosphate groups)
    phosphate_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphate_atoms:
        return False, "No phosphate group (phosphorus atom) found"

    # 2. Identify a sugar moiety by looking for a five-membered ring with 4 carbons and 1 oxygen.
    ring_info = mol.GetRingInfo()
    sugar_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
            n_oxygens = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == "O")
            n_carbons = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == "C")
            if n_oxygens == 1 and n_carbons == 4:
                sugar_ring = set(ring)
                break
    if sugar_ring is None:
        return False, "Sugar moiety (five-membered furanose ring) not found"

    # 3. Identify a nucleobase moiety:
    # Look for an aromatic ring (excluding the sugar ring) that contains at least 2 nitrogen atoms.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        # Skip if this ring is fully contained in the sugar ring.
        if sugar_ring.issuperset(ring):
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        # Check that the ring is fully aromatic.
        if not all(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        n_nitrogens = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
        if n_nitrogens >= 2:
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "Nucleobase moiety not found (no aromatic ring with sufficient nitrogen atoms)"

    # 4. Check that the phosphate group is attached to the sugar.
    # Instead of requiring that an oxygen on the phosphate is directly in the sugar ring,
    # allow for an oxygen that is directly bound to a sugar ring atom.
    phosphate_attached_to_sugar = False
    for p_atom in phosphate_atoms:
        for neighbor in p_atom.GetNeighbors():
            # Look for oxygen neighbors of P
            if neighbor.GetSymbol() != "O":
                continue
            # Direct attachment: oxygen is part of the sugar ring.
            if neighbor.GetIdx() in sugar_ring:
                phosphate_attached_to_sugar = True
                break
            # Allow attachment if the oxygen is directly bonded to a sugar ring atom.
            for o_neighbor in neighbor.GetNeighbors():
                if o_neighbor.GetIdx() in sugar_ring and o_neighbor.GetAtomicNum() in [6, 8]:
                    phosphate_attached_to_sugar = True
                    break
            if phosphate_attached_to_sugar:
                break
        if phosphate_attached_to_sugar:
            break

    if not phosphate_attached_to_sugar:
        return False, "Phosphate group is not attached to the sugar (missing condensation at the 3' or 5' OH)"
    
    # All criteria are met.
    return True, "Contains nucleoside phosphate: sugar and nucleobase moieties with phosphate attached at 3' or 5' OH"

# (Example usage - can be commented out in production)
if __name__ == "__main__":
    # Example: propanoyl-AMP 
    smiles_example = "CCC(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"
    result, reason = is_nucleotide(smiles_example)
    print(result, reason)
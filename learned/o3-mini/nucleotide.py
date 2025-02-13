"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This algorithm applies the following criteria:
  1. The molecule must contain at least one phosphorus atom.
  2. It must have a sugar moiety – we identify a five‐membered furanose ring (with 4 carbons and 1 oxygen) and then extend that set
     to include directly attached oxygen atoms (to capture the CH2OH at the 5′ position).
  3. It must contain an aromatic nucleobase (an aromatic ring with at least 2 nitrogen atoms) that is not the sugar.
  4. At least one phosphate oxygen (neighbor to P) must be “attached” to the sugar. We check that the oxygen is not part of an acyl group
     (i.e. not additionally attached to a carbonyl) and that within 1–2 bonds it reaches an atom in the sugar set.
"""

from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a nucleotide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Must contain phosphorus atom(s)
    phosphate_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphate_atoms:
        return False, "No phosphate group (P atom) found"
    
    # 2. Identify a sugar moiety.
    # We look for a 5-membered ring composed of 4 carbons and 1 oxygen.
    sugar_ring = None
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
            n_oxygens = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == "O")
            n_carbons  = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == "C")
            if n_oxygens == 1 and n_carbons == 4:
                sugar_ring = set(ring)
                break
    if sugar_ring is None:
        return False, "Sugar moiety (five‐membered furanose ring) not found"
    
    # Extend sugar group to include atoms directly attached to the furanose ring.
    # For example, the CH2OH (5′-OH) may be exocyclic.
    sugar_extended = set(sugar_ring)
    for idx in list(sugar_ring):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # We add oxygens directly attached (we could also add the corresponding carbons if desired,
            # but in our application most key attachments are via oxygen).
            if nbr.GetSymbol() == "O":
                sugar_extended.add(nbr.GetIdx())
    
    # 3. Identify an aromatic nucleobase moiety.
    # The nucleobase is defined as an aromatic ring (outside the sugar) with at least 2 nitrogen atoms.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        # Skip rings that are completely part of the sugar or sugar extension.
        if set(ring).issubset(sugar_extended):
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        n_nitrogens = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
        if n_nitrogens >= 2:
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "Nucleobase moiety not found (no aromatic ring with ≥2 nitrogen atoms outside the sugar)"
    
    # 4. Check that at least one phosphate group is attached to the sugar.
    # We iterate over each P atom, and for its oxygen neighbors we test if that O is linked to the sugar (directly,
    # or within 2 bonds) while not being part of an acyl-type connection (i.e. already bound to a carbonyl).
    phosphate_attached_to_sugar = False
    for p_atom in phosphate_atoms:
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetSymbol() != "O":
                continue
            # Exclude oxygen that is doubly bonded to a carbonyl carbon (acyl branch).
            is_acyl = False
            for nbr in o_atom.GetNeighbors():
                if nbr.GetIdx() == p_atom.GetIdx():
                    continue
                # Check if neighbor is carbon with a double bond to an oxygen.
                if nbr.GetSymbol() == "C":
                    for bond in nbr.GetBonds():
                        # If the bond to an oxygen is a double bond, mark as acyl.
                        other = bond.GetOtherAtom(nbr)
                        if other.GetSymbol() == "O" and bond.GetBondTypeAsDouble() == 2:
                            is_acyl = True
                            break
                    if is_acyl:
                        break
            if is_acyl:
                continue
            
            # Check if this oxygen (or atoms one/two bonds away) is connected to the sugar.
            # First, if the oxygen itself is in the sugar set.
            if o_atom.GetIdx() in sugar_extended:
                phosphate_attached_to_sugar = True
                break
            # Next, if any neighbor (other than the P atom) is in the sugar.
            for nbr in o_atom.GetNeighbors():
                if nbr.GetIdx() == p_atom.GetIdx():
                    continue
                if nbr.GetIdx() in sugar_extended:
                    phosphate_attached_to_sugar = True
                    break
                # Allow one more bond away: if a neighbor of this neighbor is in the sugar.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() in sugar_extended:
                        phosphate_attached_to_sugar = True
                        break
                if phosphate_attached_to_sugar:
                    break
            if phosphate_attached_to_sugar:
                break
        if phosphate_attached_to_sugar:
            break

    if not phosphate_attached_to_sugar:
        return False, "Phosphate group is not attached to the sugar (missing condensation at the 3' or 5' OH)"
    
    return True, "Contains nucleoside phosphate: sugar and nucleobase moieties with phosphate attached at 3' or 5' OH"


# (Example usage - can be commented out in production)
if __name__ == "__main__":
    # Example: propanoyl-AMP (expected to be classified as a nucleotide)
    smiles_example = "CCC(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"
    result, reason = is_nucleotide(smiles_example)
    print(result, reason)
"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This program uses several structural criteria:
1. The molecule must contain at least one phosphorus atom (P).
2. It must contain a sugar moiety: here we look for a five-membered ring (furanose) with four carbons and one oxygen.
3. It must contain a nucleobase moiety: an aromatic ring (not the sugar ring) with at least two nitrogen atoms.
4. A phosphate (P) must be attached to an oxygen that belongs to the sugar ring.
"""

from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    The algorithm checks for:
      - Presence of phosphorus (P) for the phosphate group.
      - A five-membered sugar ring (furanose) defined as a ring with 5 atoms, having 1 oxygen and 4 carbons.
      - An aromatic nucleobase ring (with at least two nitrogen atoms) that is not part of the sugar.
      - Attachment of a phosphate group to an oxygen in the sugar (indicative of condensation via the 3' or 5' OH).
    
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
    
    # 1. Check for phosphorus atoms (the presence of P suggests a phosphate group)
    phosphate_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphate_atoms:
        return False, "No phosphate group (phosphorus atom) found"
    
    # 2. Look for a sugar moiety:
    # We assume a sugar is a five-membered ring with 4 carbons and 1 oxygen.
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
    
    # 3. Look for a nucleobase moiety:
    # We require an aromatic ring (which is not the sugar) containing at least 2 nitrogen atoms.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        # Exclude rings that are the sugar ring.
        if sugar_ring.issuperset(ring):
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        # Check if every atom in the ring is aromatic
        if not all(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        # Count nitrogen atoms in this ring.
        n_nitrogens = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
        if n_nitrogens >= 2:
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "Nucleobase moiety not found (no aromatic ring with sufficient nitrogen atoms)"
    
    # 4. Check that a phosphate group is attached to the sugar:
    # Iterate over all phosphorus atoms and look for a neighboring oxygen that belongs to the sugar ring.
    phosphate_attached_to_sugar = False
    for p_atom in phosphate_atoms:
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetSymbol() == "O" and neighbor.GetIdx() in sugar_ring:
                phosphate_attached_to_sugar = True
                break
        if phosphate_attached_to_sugar:
            break
    if not phosphate_attached_to_sugar:
        return False, "Phosphate group is not attached to the sugar (missing condensation at the 3' or 5' OH)"
    
    # If all criteria are met, we call the molecule a nucleotide.
    return True, "Contains nucleoside phosphate: sugar and nucleobase moieties with phosphate attached at 3' or 5' OH"
    
# (Example usage - you may comment this out in production)
if __name__ == "__main__":
    # Test an example nucleotide: propanoyl-AMP (from user provided examples)
    smiles_example = "CCC(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"
    result, reason = is_nucleotide(smiles_example)
    print(result, reason)
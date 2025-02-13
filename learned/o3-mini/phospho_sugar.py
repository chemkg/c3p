"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: Phospho sugars 
Definition: Any monosaccharide (cyclic OR open‐chain sugar-like polyol) that contains 
an alcoholic –OH esterified with phosphoric acid.
The algorithm first looks for a phosphate group. It then checks if the oxygen that 
links the phosphate to the rest of the molecule is attached to a carbon that either 
(a) is part of a 5- or 6-membered ring containing exactly one oxygen (a typical cyclic sugar)
or 
(b) is part of a contiguous carbon chain in which at least two carbons bear hydroxyl groups.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    We require that a phosphate group is found (P with a P=O) and that one of its ester 
    oxygen atoms is attached to a carbon that is part of a sugar ring (5- or 6-membered ring with one oxygen)
    or is a part of an open-chain polyol (determined from a short DFS that finds at least 2 extra OH groups).
    
    Args:
       smiles (str): SMILES string of the molecule.
    Returns:
       bool: True if classified as phospho sugar, False otherwise.
       str: A reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- STEP 1: Identify candidate phosphate groups.
    # We iterate over phosphorus atoms that have at least one double-bonded oxygen (P=O).
    phosphate_candidates = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        has_dblO = False
        for bond in atom.GetBonds():
            # Check bond: if it is a double bond and the other atom is oxygen.
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                has_dblO = True
                break
        if has_dblO:
            phosphate_candidates.append(atom)
    
    if not phosphate_candidates:
        return False, "No phosphate group (with P=O) found in the molecule"
        
    # --- Helper function: Check if a candidate carbon is part of a sugar ring.
    def in_sugar_ring(carbon_idx):
        # Look through the rings of the molecule.
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) in (5, 6) and carbon_idx in ring:
                # Count the number of oxygen atoms in the ring.
                o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                if o_count == 1:
                    return True
        return False

    # --- Helper function: Check if a contiguous carbon-chain (open-chain sugar-like region) exists
    # starting from candidate carbon. We use a DFS of limited depth (max_depth) and count how many carbons
    # in this region (including the start) have an attached hydroxyl group (O connected with at least one hydrogen).
    def has_polyol_environment(start_atom, max_depth=3):
        visited = set()
        stack = [(start_atom, 0)]
        oh_count = 0
        
        while stack:
            atom, depth = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            # For each carbon, look for attached O bearing at least one hydrogen.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # We check explicit hydrogens: if the oxygen has any hydrogen neighbor.
                    for subnbr in nbr.GetNeighbors():
                        # Note: hydrogen atoms have atomic number 1.
                        if subnbr.GetAtomicNum() == 1:
                            oh_count += 1
                            break
            if depth < max_depth:
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                        stack.append((nbr, depth + 1))
        # A minimal polyol environment for a sugar: at least 2 OH groups (besides the one esterified).
        return oh_count >= 2

    # --- STEP 2: For each phosphate candidate, check its attached oxygens.
    for p_atom in phosphate_candidates:
        for bond in p_atom.GetBonds():
            # We skip bonds that are double bonds (P=O) so only consider ester oxygens.
            if bond.GetBondTypeAsDouble() == 2.0:
                continue
            ester_oxygen = bond.GetOtherAtom(p_atom)
            if ester_oxygen.GetAtomicNum() != 8:
                continue
            # For this ester oxygen, get the neighbor atom (other than the phosphorus).
            neighbor_atoms = [nbr for nbr in ester_oxygen.GetNeighbors() if nbr.GetIdx() != p_atom.GetIdx()]
            if not neighbor_atoms:
                continue
            # We are expecting the oxygen to be attached to a carbon
            for n_at in neighbor_atoms:
                if n_at.GetAtomicNum() != 6:
                    continue  # not a carbon, skip
                sugar_candidate = n_at
                # Option (a): if this carbon is part of a typical sugar ring.
                if in_sugar_ring(sugar_candidate.GetIdx()):
                    return True, ("Molecule contains a sugar ring (5-/6-membered with one oxygen) and a phosphate ester " 
                                  "group attached via an alcoholic oxygen")
                # Option (b): if not cyclic, check for an open-chain sugar-like (polyol) environment.
                if has_polyol_environment(sugar_candidate):
                    return True, ("Molecule contains a polyol environment (open-chain sugar-like) and a phosphate ester group " 
                                  "attached via an alcoholic oxygen")
    # If no candidate passed the tests:
    return False, "Phosphate group found but not connected to a detected sugar-like moiety"

# (Optional) For testing, one can run:
if __name__ == "__main__":
    # Feel free to test with one of the provided examples:
    test_examples = [
        ("Nc1nc(=O)[nH]c2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1", "2-hydroxy-dAMP"),
        ("O[C@H](COP(O)(O)=O)[C@@H](O)CC=O", "2-deoxy-D-ribose 5-phosphate"),
        ("O[C@H](COP(O)(O)=O)[C@@H](O)[C@H](O)[C@H](O)COP(O)(O)=O", "D-glucose 1,6-bisphosphate")
    ]
    for smi, name in test_examples:
        result, reason = is_phospho_sugar(smi)
        print(f"Example: {name}\n SMILES: {smi}\n Result: {result}\n Reason: {reason}\n")
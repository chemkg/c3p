"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: Phospho sugars 
Definition: Any monosaccharide (cyclic or open‐chain sugar‐like polyol) that contains 
an alcoholic –OH group esterified with phosphoric acid.
The algorithm first looks for a phosphate group (a phosphorus atom with at least one P=O bond). 
Then for each ester oxygen attached to that phosphorus, it checks the neighboring carbon atom.
This carbon is then examined in two ways:
  (a) If it is part of a 5- or 6-membered ring that contains exactly one ring oxygen 
      and at least two free hydroxyl groups on ring carbons and the ring has a sugar‐sized number of carbons,
  (b) Or if it is part of an open-chain (non‐cyclic) substructure made up of 3–7 carbons that carries at least two free –OH groups.
If either holds we flag the molecule as a phospho sugar.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    We require that a phosphate group exists (P with one P=O) and that one of its ester 
    oxygen atoms (not in a P=O bond) is linked to a carbon that is part of a typical sugar 
    substructure (cyclic 5-/6-membered ring or an open‐chain polyol fragment of 3–7 carbons with at least 2 free –OH groups).
    
    Args:
       smiles (str): SMILES string of the molecule.
    Returns:
       bool: True if classified as a phospho sugar, False otherwise.
       str: A reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper: returns True if an oxygen atom is a free hydroxyl group.
    # A free hydroxyl must have at least one hydrogen and should not be attached to any phosphorus.
    def is_free_hydroxyl(oxygen):
        if oxygen.GetAtomicNum() != 8:
            return False
        has_H = any(nbr.GetAtomicNum() == 1 for nbr in oxygen.GetNeighbors())
        attached_to_P = any(nbr.GetAtomicNum() == 15 for nbr in oxygen.GetNeighbors())
        return has_H and (not attached_to_P)
    
    # Helper: check if candidate carbon is part of a cyclic sugar fragment.
    # We look for rings of size 5 or 6 that contain exactly one oxygen (the ring heteroatom)
    # and (ideally) the remainder of ring atoms are carbons. We then count free –OH groups 
    # on the ring carbons (neighbors outside the ring) – requiring at least 2 free –OH groups.
    def is_cyclic_sugar(candidate_idx):
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if candidate_idx not in ring:
                continue
            if len(ring) not in (5, 6):
                continue
            # Count ring oxygens (by atomic number):
            ring_oxygens = [a_idx for a_idx in ring if mol.GetAtomWithIdx(a_idx).GetAtomicNum() == 8]
            # For a typical sugar ring (furanose or pyranose) we expect exactly one ring oxygen.
            if len(ring_oxygens) != 1:
                continue
            # Count free –OH groups on ring carbons (neighbors that are oxygen with an explicit hydrogen)
            oh_count = 0
            carbon_count = 0
            for a_idx in ring:
                atom = mol.GetAtomWithIdx(a_idx)
                if atom.GetAtomicNum() == 6:
                    carbon_count += 1
                    # Look for oxygen neighbors that are not in the ring (or not the ester oxygen)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                            if is_free_hydroxyl(nbr):
                                oh_count += 1
            # For a monosaccharide ring we expect a sugar-sized set of carbons (4 for a 5-membered ring, 5 for a 6-membered ring)
            if len(ring) == 5 and carbon_count != 4:
                continue
            if len(ring) == 6 and carbon_count != 5:
                continue
            if oh_count >= 2:
                return True
        return False

    # Helper: check if candidate carbon is part of an open-chain sugar (polyol) fragment.
    # We perform a DFS over carbon–carbon bonds and gather a connected set of carbon atoms.
    # If the group size is in the range 3 to 7 and the group has at least 2 free hydroxyl groups attached, we assume it is sugar–like.
    def is_open_chain_sugar(candidate_atom, max_depth=3):
        visited = set()
        to_visit = [(candidate_atom, 0)]
        carbon_atoms = set()
        while to_visit:
            atom, depth = to_visit.pop()
            idx = atom.GetIdx()
            if idx in visited:
                continue
            visited.add(idx)
            if atom.GetAtomicNum() == 6:
                carbon_atoms.add(idx)
                if depth < max_depth:
                    for nbr in atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                            to_visit.append((nbr, depth + 1))
        # Check if the number of carbons fits within typical monosaccharide size.
        if not (3 <= len(carbon_atoms) <= 7):
            return False
        # Count free –OH groups on all carbons in the set (neighbors attached that are oxygen and free)
        oh_count = 0
        for c_idx in carbon_atoms:
            atom = mol.GetAtomWithIdx(c_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and is_free_hydroxyl(nbr):
                    oh_count += 1
        return oh_count >= 2

    # Step 1: Identify candidate phosphate groups.
    # We iterate over phosphorus atoms that have at least one double bond (P=O) to oxygen.
    phosphate_candidates = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        has_dblO = False
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            # bond type compare: use GetBondType() and check if it equals DOUBLE
            if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_dblO = True
                break
        if has_dblO:
            phosphate_candidates.append(atom)
    
    if not phosphate_candidates:
        return False, "No phosphate group (with P=O) found in the molecule"
    
    # Step 2: For each phosphate candidate, check its attached oxygens.
    for p_atom in phosphate_candidates:
        for bond in p_atom.GetBonds():
            # Skip P=O bonds (double bonds) so only consider ester oxygens
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                continue
            ester_oxygen = bond.GetOtherAtom(p_atom)
            if ester_oxygen.GetAtomicNum() != 8:
                continue
            # For this ester oxygen, get the neighbor atoms (other than the phosphorus).
            neighbor_atoms = [nbr for nbr in ester_oxygen.GetNeighbors() if nbr.GetIdx() != p_atom.GetIdx()]
            if not neighbor_atoms:
                continue
            for nbr in neighbor_atoms:
                if nbr.GetAtomicNum() != 6:
                    continue  # we expect attachment to a carbon
                candidate_carbon = nbr
                # Check option (a): cyclic sugar substructure.
                if is_cyclic_sugar(candidate_carbon.GetIdx()):
                    return (True, "Molecule contains a sugar ring (5-/6-membered ring with one oxygen and at least 2 free OH groups) and a phosphate ester group attached via an alcoholic oxygen")
                # Check option (b): open-chain sugar-like (polyol) environment.
                if is_open_chain_sugar(candidate_carbon):
                    return (True, "Molecule contains an open-chain sugar-like (polyol) environment (3–7 carbons with at least 2 free OH groups) and a phosphate ester group attached via an alcoholic oxygen")
    
    # If no candidate passed the tests:
    return False, "Phosphate group found but not connected to a detected monosaccharide-like moiety"

# (Optional) For testing, one can run:
if __name__ == "__main__":
    # Test with a few examples given in the assignment
    test_examples = [
        ("Nc1nc(=O)[nH]c2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1", "2-hydroxy-dAMP"),
        ("O[C@H](COP(O)(O)=O)[C@@H](O)CC=O", "2-deoxy-D-ribose 5-phosphate"),
        ("OCC(=O)[C@H](O)[C@H](O)[C@@H](O)COP(O)(O)=O", "keto-L-tagatose 6-phosphate"),
        ("[H]C(=O)[C@H](O)[C@H](O)COP(O)(O)=O", "D-erythrose 4-phosphate")
    ]
    for smi, name in test_examples:
        result, reason = is_phospho_sugar(smi)
        print(f"Example: {name}\n SMILES: {smi}\n Result: {result}\n Reason: {reason}\n")
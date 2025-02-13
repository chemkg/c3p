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
      and at least one free hydroxyl group on one of the ring carbons (other than the phosphorylated position)
      and the ring has a sugar‐sized number of carbons (4 for furanoses, 5 for pyranoses),
  (b) Or if it is part of an open‐chain (non‐cyclic) fragment made up of 3–7 carbons (only non‐ring carbons)
      that carries at least one free –OH group.
If either holds we flag the molecule as a phospho sugar.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar must contain a phosphate group (P with at least one P=O) 
    and one of its ester oxygen atoms (i.e. an O attached via a single bond to phosphorus) 
    must be connected to a carbon that is part of a typical sugar (cyclic 5-/6-membered ring with one ring oxygen 
    and at least one free OH group; or an open‐chain polyol fragment of 3–7 non‐ring carbons having at least one free OH group).
    
    Args:
       smiles (str): SMILES string of the molecule.
    Returns:
       bool: True if the molecule is classified as a phospho sugar, False otherwise.
       str: A reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper: returns True if an oxygen atom is considered a free hydroxyl.
    # We require that the oxygen has at least one hydrogen and is not attached to any phosphorus.
    def is_free_hydroxyl(oxygen):
        if oxygen.GetAtomicNum() != 8:
            return False
        # Check that at least one neighbor is H
        has_H = any(nbr.GetAtomicNum() == 1 for nbr in oxygen.GetNeighbors())
        # Ensure not directly attached to phosphorus (aside from the ester bond under consideration)
        attached_to_P = any(nbr.GetAtomicNum() == 15 for nbr in oxygen.GetNeighbors())
        return has_H and (not attached_to_P)

    # Helper: Determine if a candidate carbon is in a cyclic sugar fragment.
    # We expect a ring of 5 or 6 atoms that contains exactly one oxygen atom.
    # The sugar ring normally has one less carbon than the ring size; for a 5-membered ring (furanose) expect 4 carbons; 
    # for 6-membered ring (pyranose) expect 5 carbons.
    # We also now require at least one free hydroxyl (other than the one bearing the phosphate) on any carbon of the ring.
    def is_cyclic_sugar(candidate_idx):
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if candidate_idx not in ring:
                continue
            if len(ring) not in (5, 6):
                continue
            # Count ring oxygens (atomic number 8):
            ring_oxygens = [a_idx for a_idx in ring if mol.GetAtomWithIdx(a_idx).GetAtomicNum() == 8]
            if len(ring_oxygens) != 1:
                continue
            # Count how many ring carbons are present (should be 4 if ring size 5, and 5 if ring size 6)
            ring_carbons = [a_idx for a_idx in ring if mol.GetAtomWithIdx(a_idx).GetAtomicNum() == 6]
            if len(ring) == 5 and len(ring_carbons) != 4:
                continue
            if len(ring) == 6 and len(ring_carbons) != 5:
                continue
            # Count free hydroxyl groups on ring carbons. Skip the carbon already known to be connected to phosphate.
            oh_count = 0
            for a_idx in ring:
                atom = mol.GetAtomWithIdx(a_idx)
                if atom.GetAtomicNum() != 6:
                    continue
                # Check neighbors that are oxygen and are not part of the ring.
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                        if is_free_hydroxyl(nbr):
                            oh_count += 1
            if oh_count >= 1:
                return True
        return False

    # Helper: Determine if a candidate carbon is part of an open-chain sugar-like fragment.
    # We perform a DFS over connected non-ring carbons.
    def is_open_chain_sugar(candidate_atom):
        # Only consider candidate atom if it is not in any ring.
        if mol.GetAtomWithIdx(candidate_atom.GetIdx()).IsInRing():
            return False
        visited = set()
        stack = [candidate_atom]
        carbon_ids = set()
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6 and (not atom.IsInRing()):
                carbon_ids.add(atom.GetIdx())
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()) and nbr.GetIdx() not in visited:
                        stack.append(nbr)
        # Check that the number of contiguous carbons is in a monosaccharide-sized range
        if not (3 <= len(carbon_ids) <= 7):
            return False
        # Count free hydroxyl groups on these carbons (neighbors that are oxygen and free)
        oh_count = 0
        for cid in carbon_ids:
            atom = mol.GetAtomWithIdx(cid)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and is_free_hydroxyl(nbr):
                    oh_count += 1
        return oh_count >= 1

    # Step 1: Identify candidate phosphate groups.
    # Look for phosphorus atoms with at least one double bond to an oxygen.
    phosphate_candidates = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        has_dblO = False
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_dblO = True
                break
        if has_dblO:
            phosphate_candidates.append(atom)
    
    if not phosphate_candidates:
        return False, "No phosphate group (with at least one P=O bond) found in the molecule"
    
    # Step 2: For each candidate phosphate, examine its ester oxygen atoms.
    for p_atom in phosphate_candidates:
        for bond in p_atom.GetBonds():
            # Skip double bonds (P=O) so we only consider ester bonds (single bonds)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                continue
            ester_oxygen = bond.GetOtherAtom(p_atom)
            if ester_oxygen.GetAtomicNum() != 8:
                continue
            # Find neighbor atoms of the oxygen (other than phosphorus)
            neighbor_atoms = [nbr for nbr in ester_oxygen.GetNeighbors() if nbr.GetIdx() != p_atom.GetIdx()]
            if not neighbor_atoms:
                continue
            for nbr in neighbor_atoms:
                if nbr.GetAtomicNum() != 6:
                    continue  # We expect the phosphate to be attached via an oxygen to a carbon.
                candidate_carbon = nbr
                # Option (a): Check if the carbon is part of a cyclic sugar fragment.
                if is_cyclic_sugar(candidate_carbon.GetIdx()):
                    return True, ("Molecule contains a sugar ring (5- or 6-membered ring with one oxygen and at least one free OH group) " 
                                "and a phosphate ester group attached via an alcoholic oxygen")
                # Option (b): Check if the carbon is part of an open-chain sugar fragment.
                if is_open_chain_sugar(candidate_carbon):
                    return True, ("Molecule contains an open-chain sugar-like fragment (3–7 non-ring carbons with at least one free OH group) " 
                                 "and a phosphate ester group attached via an alcoholic oxygen")
    
    return False, "Phosphate group found but not connected to a detected monosaccharide-like moiety"

# Optional test block:
if __name__ == "__main__":
    # List of test SMILES from the assignment with names for clarity
    test_examples = [
        ("Nc1nc(=O)[nH]c2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1", "2-hydroxy-dAMP"),
        ("O[C@H](COP(O)(O)=O)[C@@H](O)CC=O", "2-deoxy-D-ribose 5-phosphate"),
        ("CN(C)c1nc(N)c2ncn([C@@H]3O[C@@H]4COP(O)(=O)O[C@H]4[C@H]3O)c2n1", "2-(dimethylamino)-cAMP"),
        ("OCC(=O)[C@H](O)[C@H](O)[C@@H](O)COP(O)(O)=O", "keto-L-tagatose 6-phosphate"),
        ("[H]C(=O)[C@H](O)[C@H](O)COP(O)(O)=O", "D-erythrose 4-phosphate")
    ]
    for smi, name in test_examples:
        result, reason = is_phospho_sugar(smi)
        print(f"Example: {name}\n SMILES: {smi}\n Result: {result}\n Reason: {reason}\n")
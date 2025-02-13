"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: Phytosterols
Definition: Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.
This program attempts to detect a steroid nucleus (gonane: 3 fused six‐membered rings and 1 fused five‐membered ring with 17 carbons),
a hydroxyl group (either free –OH or present as a glycosidic substituent) attached to a nucleus carbon,
and a side chain (a carbon branch from the nucleus that extends for a few carbons).
Note: This is a heuristic approach and may produce false positives/negatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    The program performs the following:
      (1) Looks for the tetracyclic steroid (gonane) nucleus as a fused system of 3 six‐membered rings and 1 five‐membered ring with 17 atoms.
      (2) Checks for a free –OH group attached directly to a nucleus carbon. If not found, it attempts to detect an oxygen
          attached to the nucleus that may be part of a glycoside.
      (3) Checks for a side chain attached to one of the nucleus atoms that has at least a minimal (here, 3–4 carbon) chain.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        (bool, str): True with explanation if the molecule is classified as a phytosterol,
                     otherwise False with reason for rejection.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information to look for fused rings.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Identify rings of size 6 and 5.
    six_rings = [r for r in ring_info if len(r) == 6]
    five_rings = [r for r in ring_info if len(r) == 5]
    
    nucleus_atom_indices = None
    # Look for a combination of 3 six‐membered rings and 1 five‐membered ring whose union gives 17 atoms
    for six_combo in itertools.combinations(six_rings, 3):
        for five_ring in five_rings:
            union_atoms = set()
            for ring in six_combo:
                union_atoms.update(ring)
            union_atoms.update(five_ring)
            if len(union_atoms) == 17:
                nucleus_atom_indices = union_atoms
                break
        if nucleus_atom_indices is not None:
            break
    
    if nucleus_atom_indices is None:
        return False, "Steroid (gonane) nucleus not found (requires 3 six‐membered and 1 five‐membered fused rings with 17 atoms)"
    
    # ---------- Check for a hydroxyl group attached to the nucleus ----------
    # First attempt: look for free hydroxyl groups using the SMARTS [OX2H]
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    hydroxyl_in_nucleus = False
    for match in oh_matches:
        oh_atom_idx = match[0]
        oh_atom = mol.GetAtomWithIdx(oh_atom_idx)
        # Look if any neighbor is a carbon in the nucleus.
        for neigh in oh_atom.GetNeighbors():
            if neigh.GetAtomicNum() == 6 and neigh.GetIdx() in nucleus_atom_indices:
                hydroxyl_in_nucleus = True
                break
        if hydroxyl_in_nucleus:
            break

    # Second attempt: if no free hydroxyl was found, check for an oxygen attached to the nucleus that could be part of a glycoside.
    # Heuristic: if an oxygen is attached to a nucleus carbon and that oxygen itself is bonded to at least one other oxygen,
    # we consider that as a potential “masked” hydroxyl.
    if not hydroxyl_in_nucleus:
        for idx in nucleus_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6: 
                continue
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in nucleus_atom_indices:
                    # Check if this O neighbor has other oxygen neighbors (a common pattern in sugars)
                    o_neighbors = [n for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 8]
                    if len(o_neighbors) >= 1:
                        hydroxyl_in_nucleus = True
                        break
            if hydroxyl_in_nucleus:
                break

    if not hydroxyl_in_nucleus:
        return False, "No hydroxyl (–OH) group (free or as part of a glycoside) found attached to the steroid nucleus"
    
    # ---------- Check for a side chain attached to the nucleus ----------
    # We look at each nucleus atom and see if it is attached to a carbon not part of the nucleus.
    # Then we try to follow a carbon chain (using DFS) to see if it extends for at least 3 additional carbons.
    def chain_length(atom_idx, visited):
        """Recursively count the length of a consecutive carbon chain starting from atom_idx.
           Only traverse carbon atoms not in nucleus and not visited before.
        """
        mol_atom = mol.GetAtomWithIdx(atom_idx)
        max_length = 0
        for nbr in mol_atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in visited:
                continue
            if nbr.GetAtomicNum() == 6 and n_idx not in nucleus_atom_indices:
                visited.add(n_idx)
                branch_length = 1 + chain_length(n_idx, visited)
                visited.remove(n_idx)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    side_chain_found = False
    for idx in nucleus_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in nucleus_atom_indices:
                clength = chain_length(nbr.GetIdx(), set())
                # Here we require at least 3 carbons in the side chain (which is a weak but reasonable rule).
                if clength >= 3:
                    side_chain_found = True
                    break
        if side_chain_found:
            break
            
    if not side_chain_found:
        return False, "No adequate carbon side chain detected attached to the steroid nucleus"
    
    return True, "Molecule contains a steroid nucleus (gonane) with an –OH (free or masked) and an extended side chain, consistent with phytosterols"

# Example usage:
if __name__ == "__main__":
    # Test a couple of examples (one known phytosterol and one borderline case)
    examples = {
        "ergosterol": "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)",
        "stigmasterol 3-O-beta-D-glucoside": "[H][C@@]1(O[C@H]2CC[C@@]3(C)C(C2)=CC[C@]2([H])[C@]3([H])CC[C@]3(C)[C@](C)(CC[C@@H](C)C(C)C)O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)"
    }
    for name, smi in examples.items():
        res, reason = is_phytosterols(smi)
        print(f"{name}: {res}, {reason}")
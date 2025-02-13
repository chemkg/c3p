"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: Phytosterols
Definition: Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.
This heuristic program detects a steroid nucleus (gonane) by looking for a set of three fused six‐membered rings plus one five‐membered ring (with 17 unique atoms),
then checks that at least one oxygen substituent attached to a nucleus carbon appears as a free hydroxyl or a hydroxyl that is “masked” (for example, in a glycoside).
It also requires that the nucleus carries a carbon side chain of minimum length.
Some non‐phytosterol steroids may be misclassified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    The algorithm attempts to detect:
     (1) A steroid (gonane) nucleus as a fused system of three six‐membered rings and one five‐membered ring (union of 17 atoms)
     (2) A hydroxyl function attached (directly or through a glycosidic substituent) to the nucleus.
         (It avoids groups that appear to be acylated (e.g. acetoxy)).
     (3) A side chain attached to a nucleus carbon that extends for at least 3 carbon atoms.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        tuple(bool, str): True with explanation if classified as phytosterol,
                          otherwise False with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- STEP 1: Identify the steroid nucleus (gonane) ---
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Partition rings by size; we are interested in 6-membered and 5-membered rings
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
    
    # --- STEP 2: Check for a hydroxyl substituent on the nucleus ---
    # We want to see a free –OH or one that is part of a glycoside attached to a nucleus carbon.
    # Also we want to avoid acyl (ester) groups.
    
    # Helper: detect if an oxygen atom is part of an acyl group (i.e. bonded to a carbon that is doubly bonded to another oxygen)
    def is_acyl(oxygen_atom):
        for nbr in oxygen_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                for nn in nbr.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                    if nn.GetAtomicNum() == 8 and bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        return True
        return False
    
    # Helper: given a starting atom (outside the nucleus), return the connected fragment (set of atom indices) not including ones in the nucleus.
    def get_substituent_fragment(start_idx, nucleus_set):
        visited = set()
        stack = [start_idx]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            for nbr in atom.GetNeighbors():
                n_idx = nbr.GetIdx()
                if n_idx not in nucleus_set and n_idx not in visited:
                    stack.append(n_idx)
        return visited

    hydroxyl_found = False
    # For each nucleus atom, look at neighbors (substituents) that are not in the nucleus.
    for nuc_idx in nucleus_atom_indices:
        nuc_atom = mol.GetAtomWithIdx(nuc_idx)
        for nbr in nuc_atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in nucleus_atom_indices:
                continue
            # Check target atom if it is oxygen
            if nbr.GetAtomicNum() == 8:
                # If this oxygen is part of an acyl (ester) group, skip it.
                if is_acyl(nbr):
                    continue
                # If this oxygen has an explicit hydrogen, count as hydroxyl.
                if nbr.GetTotalNumHs() > 0:
                    hydroxyl_found = True
                    break
                else:
                    # Otherwise, check the substituent fragment attached from this oxygen
                    frag = get_substituent_fragment(n_idx, nucleus_atom_indices)
                    # Look for any oxygen in that fragment that carries at least one hydrogen.
                    for idx in frag:
                        atom_frag = mol.GetAtomWithIdx(idx)
                        if atom_frag.GetAtomicNum() == 8 and atom_frag.GetTotalNumHs() > 0:
                            hydroxyl_found = True
                            break
                if hydroxyl_found:
                    break
        if hydroxyl_found:
            break

    if not hydroxyl_found:
        return False, "No hydroxyl (–OH) group (free or as part of a glycoside) found attached to the steroid nucleus"
    
    # --- STEP 3: Check for the presence of an extended carbon side chain ---
    # We search for a side chain attached to any nucleus atom that is a carbon not in the nucleus.
    # We then recursively follow chains (only carbon atoms) and require at least 3 carbons in length.
    def chain_length(atom_idx, visited):
        """Recursively count the length of a carbon-only chain starting from atom_idx.
           Do not traverse back into the nucleus.
        """
        atom = mol.GetAtomWithIdx(atom_idx)
        max_len = 0
        for nbr in atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in visited or n_idx in nucleus_atom_indices:
                continue
            if nbr.GetAtomicNum() == 6:
                visited.add(n_idx)
                branch = 1 + chain_length(n_idx, visited)
                visited.remove(n_idx)
                if branch > max_len:
                    max_len = branch
        return max_len

    side_chain_found = False
    for nuc_idx in nucleus_atom_indices:
        nuc_atom = mol.GetAtomWithIdx(nuc_idx)
        for nbr in nuc_atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in nucleus_atom_indices:
                continue
            if nbr.GetAtomicNum() == 6:  # potential start of side chain
                clength = chain_length(n_idx, set([nuc_idx]))
                if clength >= 3:
                    side_chain_found = True
                    break
        if side_chain_found:
            break

    if not side_chain_found:
        return False, "No adequate carbon side chain detected attached to the steroid nucleus"
    
    # --- If all criteria are met, we classify as phytosterol ---
    return True, ("Molecule contains a steroid nucleus (gonane) with an –OH (either free or in a glycoside) "
                  "and an extended aliphatic side chain, consistent with phytosterols")

# Example usage:
if __name__ == "__main__":
    examples = {
        "ergosterol": "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)\\C=C\\[C@H](C)C(C)C",
        "stigmasterol 3-O-beta-D-glucoside": "[H][C@@]1(O[C@H]2CC[C@@]3(C)C(C2)=CC[C@]2([H])[C@]3([H])CC[C@]3(C)[C@](C)(CC[C@@H](C)C(C)C)O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)",
        "3beta-hydroxy-9beta-9,19-cyclolanost-24-en-28-oic acid": "[C@]123[C@@]4([C@]([C@]([C@@H](O)CC4)(C)C(=O)O)(CC[C@]1([C@]5(C)CC[C@@]([C@]5(CC2)C)([C@@H](CCC=C(C)C)C)[H])[H])[H])C3"
    }
    for name, smi in examples.items():
        res, reason = is_phytosterols(smi)
        print(f"{name}: {res}, {reason}")
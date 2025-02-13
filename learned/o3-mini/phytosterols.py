"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: Phytosterols
Definition: Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.
This program checks for the steroid nucleus (the “gonane” system of three cyclohexanes and one cyclopentane fused together with 17 skeletal carbons),
the presence of a hydroxyl group on the nucleus, and an attached carbon side chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    The program searches for:
      (1) A steroid nucleus: the tetracyclic gonane framework consisting of 3 six‐membered rings and 1 five‐membered ring
          that share a total of 17 carbon atoms.
      (2) A hydroxyl group attached to one of the nucleus atoms (typically at C3).
      (3) A side chain attached to the nucleus (e.g. at C17).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a phytosterol, False otherwise.
        str: Explanation/reason for the classification decision.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute ring information for the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # In the gonane (steroid) nucleus we expect 4 fused rings:
    # exactly 3 six-membered rings and 1 five-membered ring; the fused set should contain 17 atoms.
    six_rings = [r for r in ring_info if len(r) == 6]
    five_rings = [r for r in ring_info if len(r) == 5]
    
    nucleus_atom_indices = None
    # Try all combinations: pick 3 six-membered rings and 1 five-membered ring.
    for six_combo in itertools.combinations(six_rings, 3):
        for five_ring in five_rings:
            # Union of atom indices from these 4 rings.
            union_atoms = set()
            for ring in six_combo:
                union_atoms.update(ring)
            union_atoms.update(five_ring)
            # In the classical steroid nucleus the 4 fused rings share atoms so that the total number is 17.
            if len(union_atoms) == 17:
                nucleus_atom_indices = union_atoms
                break
        if nucleus_atom_indices is not None:
            break
    
    if nucleus_atom_indices is None:
        return False, "Steroid (gonane) nucleus not found (requires 3 six-membered and 1 five-membered fused rings with 17 atoms)"
        
    # Check for a hydroxyl group attached to an atom that is part of the nucleus.
    # We look for an -OH group using the SMARTS pattern [OX2H].
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    hydroxyl_in_nucleus = False
    for match in oh_matches:
        # match is a tuple of atom indices matching the oxygen atom.
        if match[0] in nucleus_atom_indices:
            hydroxyl_in_nucleus = True
            break
    if not hydroxyl_in_nucleus:
        return False, "No hydroxyl (–OH) group found on the steroid nucleus"
    
    # Check for a side chain: we expect at least one carbon neighbor attached to one of the nucleus atoms that is outside the nucleus.
    side_chain_found = False
    for idx in nucleus_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in nucleus_atom_indices and neighbor.GetAtomicNum() == 6:
                side_chain_found = True
                break
        if side_chain_found:
            break
    if not side_chain_found:
        return False, "No carbon side chain detected attached to the steroid nucleus"
    
    # Optionally, one can check for additional features such as overall molecular weight or
    # further details of the side chain. Here we assume that the presence of the gonane nucleus,
    # the hydroxyl group, and an external carbon chain is enough.
    
    # All tests passed – the molecule is consistent with being a phytosterol.
    return True, "Molecule contains a steroid nucleus (gonane), a hydroxyl group on the nucleus, and a side chain, consistent with phytosterols"

# Example usage (uncomment to test):
# smiles_examples = [
#     # ergosterol, campesterol, stigmasterol, etc.
#     "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)\\C=C\\[C@H](C)C(C)C",
# ]
# for s in smiles_examples:
#     res, reason = is_phytosterols(s)
#     print(res, reason)
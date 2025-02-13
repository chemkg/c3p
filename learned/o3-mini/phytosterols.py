"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: Phytosterols
Definition: Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.
This program checks for the tetracyclic steroid (gonane) nucleus, a hydroxyl group attached to a nucleus carbon,
and a carbon side chain attached to the nucleus.
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
      (2) A hydroxyl group attached to one of the nucleus carbons.
      (3) A side chain attached to the nucleus (e.g. at C17).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a phytosterol, False otherwise.
        str: Explanation/reason for classification decision.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # For steroids, we expect 3 six-membered rings and 1 five-membered ring that are fused
    # together so that the union of their atoms equals 17. This will be our steroid nucleus.
    six_rings = [r for r in ring_info if len(r) == 6]
    five_rings = [r for r in ring_info if len(r) == 5]
    
    nucleus_atom_indices = None
    # Try all combinations of 3 six-membered rings and 1 five-membered ring to find a fused set with 17 atoms.
    for six_combo in itertools.combinations(six_rings, 3):
        for five_ring in five_rings:
            union_atoms = set()
            for ring in six_combo:
                union_atoms.update(ring)
            union_atoms.update(five_ring)
            if len(union_atoms) == 17:  # Found the gonane skeleton
                nucleus_atom_indices = union_atoms
                break
        if nucleus_atom_indices is not None:
            break
    
    if nucleus_atom_indices is None:
        return False, "Steroid (gonane) nucleus not found (requires 3 six-membered and 1 five-membered fused rings with 17 atoms)"
        
    # Check for a hydroxyl group (–OH) attached to a nucleus carbon.
    # The SMARTS [OX2H] will pick up hydroxyl groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    hydroxyl_in_nucleus = False
    for match in oh_matches:
        oh_atom_index = match[0]
        oh_atom = mol.GetAtomWithIdx(oh_atom_index)
        # Check if any neighbor of the hydroxyl oxygen is a carbon that belongs to the nucleus.
        for neighbor in oh_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in nucleus_atom_indices:
                hydroxyl_in_nucleus = True
                break
        if hydroxyl_in_nucleus:
            break

    if not hydroxyl_in_nucleus:
        return False, "No hydroxyl (–OH) group found attached to the steroid nucleus"
    
    # Check for a side chain: look for a carbon neighbor attached to one of the nucleus atoms that lies outside the nucleus.
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
    
    # If all criteria are met, classify as a phytosterol.
    return True, "Molecule contains a steroid nucleus (gonane), a hydroxyl group attached to the nucleus, and a side chain, consistent with phytosterols"

# Uncomment the example usage to test:
# smiles_examples = [
#     "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)\\C=C\\[C@H](C)C(C)C",  # ergosterol
#     # Other examples can be added here...
# ]
# for s in smiles_examples:
#     res, reason = is_phytosterols(s)
#     print(res, reason)
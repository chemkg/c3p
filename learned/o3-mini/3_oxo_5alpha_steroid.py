"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: 3-oxo-5alpha-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5.
Heuristic: 
  - Look for a ketone group (C=O) that is part of a ring.
  - Check that the molecule contains four fused rings (typical steroid nucleus: three six-membered rings and one five-membered ring).
  - Among the atoms in these rings, look for a chiral center located at a ring junction (i.e. present in at least two of the steroid rings) with the stereochemistry designated as '@@' (a heuristic for an alpha oriented hydrogen).
Note: Steroid stereochemistry is complex. If the analysis is inconclusive or for molecules that do not clearly satisfy these criteria, classification may be incorrect.
"""

from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a 3-oxo-5alpha-steroid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # --- Step 1: Check for the 3-oxo function ---
    # We require a ketone [C]=O group in a ring.
    ketone_smarts = "[#6;R](=O)"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found in a ring (3-oxo requirement not met)."
    
    # --- Step 2: Check for the steroid nucleus ---
    # A typical steroid nucleus consists of four fused rings:
    # three six-membered rings and one five-membered ring.
    # We use a simple heuristic: count rings with 5 or 6 members.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples of atom indices forming rings
    steroid_rings = [ring for ring in atom_rings if len(ring) in (5, 6)]
    if len(steroid_rings) < 4:
        return False, "Fewer than four rings of size 5 or 6 found (steroid nucleus not detected)."
    
    # --- Step 3: Check for 5alpha configuration ---
    # For a 5α steroid we expect the atom at position 5 (typically a ring junction in the steroid nucleus)
    # to be chiral with an alpha configuration.
    # Here we look for a chiral center (as reported by RDKit) that is present in at least 2 of the steroid rings
    # and has stereochemistry specified as '@@' (a heuristic for alpha configuration).
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    five_alpha_found = False
    for atom_idx, chir in chiral_centers:
        # Count how many steroid rings (from our subset) include this atom.
        n_in_rings = sum(1 for ring in steroid_rings if atom_idx in ring)
        if n_in_rings >= 2:
            # Heuristic: if the chiral tag is '@@', we take it to mean alpha oriented.
            if chir == "@@":
                five_alpha_found = True
                break
    if not five_alpha_found:
        return False, "No chiral center with alpha (\"@@\") configuration at a ring junction found (5α requirement not met)."

    return True, "Molecule has a 3-oxo group in a steroid nucleus and a chiral center with alpha configuration (likely 5α)."

# Example usage (for testing; remove or comment out in production):
if __name__ == "__main__":
    # Test with one of the provided examples: 5alpha-Androstane-3,11,17-trione
    test_smiles = "C[C@]12CC(=O)[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2=O"
    result, reason = is_3_oxo_5alpha_steroid(test_smiles)
    print(result, reason)
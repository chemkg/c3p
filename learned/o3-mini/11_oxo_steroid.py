"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: CHEBI/Custom: 11-oxo steroid
Definition: An oxo steroid must have a fused steroid nucleus with exactly 17 carbons 
            and at least one ring-bound ketone group (C(=O)) situated in that nucleus.
Heuristic implementation:
  1. Parse the SMILES.
  2. Identify all rings that contain only carbons.
  3. Build a “fusion” graph where rings are connected if they share at least 2 atoms.
  4. From the connected components, select the one with the most nucleus atoms.
  5. Require that the nucleus (the union of atoms in that fused component) contains exactly 17 carbon atoms.
  6. Then, look for a candidate ketone:
       - It must be a carbon atom in the nucleus that is part of at least 2 rings from that nucleus.
       - It must be double-bonded to an oxygen.
       - Its two other neighbors (aside from the oxygen) must be carbons in the nucleus.
       - At least one of the rings containing this candidate should be 6-membered.
       
  If these conditions are met, we classify the molecule as an 11-oxo steroid.
  
Note: This heuristic does not perform full stereochemical assignment or explicit numbering.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    The heuristic:
      - Parse the molecule.
      - Identify all carbon-only rings and build a fusion graph (rings are fused if they share at least 2 atoms).
      - Determine the largest connected component (the steroid nucleus candidate) and require it
        to have exactly 17 carbon atoms.
      - Within that nucleus, search for a ketone group (C(=O)) on a carbon that belongs to at least 
        two nucleus rings, with both additional substituents being carbons of the nucleus. Also, require 
        at least one associated ring to be six-membered (since position 11 lies on a six-membered ring).
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule is classified as an 11-oxo steroid, False otherwise.
      str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry and compute 2D coordinates (helps with ring perception)
    Chem.AssignStereochemistry(mol, cleanIt=True)
    AllChem.Compute2DCoords(mol)
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    # Find rings that are exclusively carbons
    carbon_rings = []
    for ring in all_rings:
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            # Save as set for easy intersection computation and preserve the ring size
            carbon_rings.append(set(ring))
    
    if not carbon_rings:
        return False, "No carbon-only rings found; not a steroid-like nucleus."
    
    # Build a graph among rings: rings are fused if they share at least 2 atoms.
    n_rings = len(carbon_rings)
    graph = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(carbon_rings[i].intersection(carbon_rings[j])) >= 2:
                graph[i].add(j)
                graph[j].add(i)
                
    # Find connected components of the ring graph via DFS.
    seen = set()
    components = []
    for i in range(n_rings):
        if i in seen:
            continue
        stack = [i]
        comp = set()
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp.add(cur)
            for nbr in graph[cur]:
                if nbr not in comp:
                    stack.append(nbr)
        seen.update(comp)
        components.append(comp)
    
    # For each connected component, take the union of all atom indices in its rings.
    nucleus_candidates = []
    for comp in components:
        nucleus_atoms = set()
        for ring_idx in comp:
            nucleus_atoms.update(carbon_rings[ring_idx])
        # Count how many carbons in this nucleus.
        nucleus_carbon_count = sum(1 for idx in nucleus_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        nucleus_candidates.append((comp, nucleus_atoms, nucleus_carbon_count))
    
    # Choose the candidate with the most nucleus atoms.
    nucleus_candidates.sort(key=lambda x: x[2], reverse=True)
    if not nucleus_candidates:
        return False, "No fused ring nucleus candidate found."
    chosen_comp, nucleus_atom_idxs, nucleus_carbon_count = nucleus_candidates[0]
    
    # For a typical steroid nucleus, we expect exactly 17 carbon atoms.
    if nucleus_carbon_count != 17:
        return False, f"Fused nucleus has {nucleus_carbon_count} carbons; expected 17 for a steroid nucleus."
    
    # For each ring in the chosen component, record its size (for later reference)
    nucleus_rings = []
    for ring_idx in chosen_comp:
        ring_set = carbon_rings[ring_idx]
        # Only consider rings fully contained in the nucleus (they should be, by construction)
        if ring_set.issubset(nucleus_atom_idxs):
            nucleus_rings.append((ring_idx, ring_set, len(ring_set)))
    
    # Now search for a candidate ketone group within the nucleus.
    # The ketone candidate must:
    #   - be a carbon (atomic number 6) present in the nucleus
    #   - be bonded by a double bond to an oxygen (atomic number 8)
    #   - have exactly two other neighbors; both must be carbons and in the nucleus.
    #   - be in at least two nucleus rings; and at least one of these rings should be six-membered.
    ketone_found = False
    ketone_details = ""
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetIdx() not in nucleus_atom_idxs:
            continue
        # Find all nucleus rings (by index) that contain this atom.
        rings_here = []
        for ring_idx, ring_set, ring_size in nucleus_rings:
            if atom.GetIdx() in ring_set:
                rings_here.append((ring_idx, ring_size))
        if len(rings_here) < 2:
            continue  # must belong to at least 2 rings in the nucleus
        
        # Look for a double bond to oxygen.
        for bond in atom.GetBonds():
            # Check the bond is a double bond.
            if bond.GetBondTypeAsDouble() != 2:
                continue
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() != 8:
                continue
            # Now check that aside from this oxygen, the carbonyl carbon has exactly 2 other neighbors.
            other_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() != other.GetIdx()]
            if len(other_neighbors) != 2:
                continue
            # Both other neighbors must be carbons in the nucleus.
            if not all(nbr.GetAtomicNum() == 6 and nbr.GetIdx() in nucleus_atom_idxs for nbr in other_neighbors):
                continue
            # Also check that at least one of the rings (among those that contain the ketone carbon) is six-membered.
            if not any(ring_size == 6 for (_, ring_size) in rings_here):
                continue
            # If all criteria are met, we presume this is the 11-oxo.
            ketone_found = True
            ketone_details = (f"Atom {atom.GetIdx()} forms a ketone (C=O with atom {other.GetIdx()}) "
                              f"and is in {len(rings_here)} nucleus rings (ring sizes: " +
                              ", ".join(str(ring_size) for (_, ring_size) in rings_here) + ")")
            break
        if ketone_found:
            break

    if not ketone_found:
        return False, "No suitable ring-bound ketone (11-oxo) found in the steroid nucleus."
    
    reason = ("Molecule contains a steroid nucleus with 17 carbons and a candidate 11-oxo group: " +
              ketone_details)
    return True, reason

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles_list = [
        # Emesterone B
        "O=C1C=C2C=CC3=C4[C@]5([C@@H]([C@H](C)[C@@H](C5)[C@@H](O)[C@@H](C(C)C)C)CC4)CC([C@]3([C@@]2(C)CC1)O)=O",
        # 6beta-hydroxycortisone
        "C1[C@@]2([C@@]([C@@]3(C([C@@H]1O)=CC(CC3)=O)C)(C(C[C@@]4([C@@](CC[C@@]24[H])(O)C(CO)=O)C)=O)[H])[H]"
    ]
    for smi in test_smiles_list:
        result, detail = is_11_oxo_steroid(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Detail:", detail)
        print("--------------------------------------------------")
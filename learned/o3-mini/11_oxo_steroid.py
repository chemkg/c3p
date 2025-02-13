"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: CHEBI/Custom: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11.
Heuristic implementation:
  1. The SMILES is parsed.
  2. We require that the molecule contains a fused steroid nucleus – here defined
     as a set of four fused “carbon‐only” rings (typically one five-membered ring and three six-membered rings,
     or four six‐membered rings) that share bonds (we require an overlap of at least two atoms between rings).
  3. The molecule must contain at least one candidate ring‐bound ketone group.
     For a ketone (C(=O)) to be acceptable:
       (a) the carbonyl carbon must be in the fused nucleus (i.e. belong to at least 2 rings in the nucleus),
       (b) the bond to oxygen must be a double bond (forming a C=O),
       (c) the two other substituents on that carbon are carbons and are in the fused nucleus.
  4. If these criteria are met, then we classify the molecule as an 11-oxo steroid.
Note: This heuristic does not perform full stereochemical or explicit numbering assignment.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    The heuristic:
      - Parse the molecule.
      - Determine all carbon-only rings and find connected sets of rings that are fused
        (sharing at least 2 atoms). Accept only if one connected component has exactly 4 rings.
      - Identify at least one candidate ketone group (C(=O)) on a carbon that is in at least 2 rings
        of the fused nucleus and whose other two substituents are carbons in the fused nucleus.
    
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
    
    # Ensure stereochemistry assignment and ring perception
    Chem.AssignStereochemistry(mol, cleanIt=True)
    AllChem.Compute2DCoords(mol)
    
    # Gather ring information and restrict to rings with only carbon atoms.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    
    carbon_rings = []
    for ring in all_rings:
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            carbon_rings.append(set(ring))  # store as set for intersection tests
            
    if not carbon_rings:
        return False, "No carbon-only rings found; not a steroid-like nucleus."
    
    # Build connectivity graph among rings: two rings are fused if they share at least 2 atoms.
    # Represent the graph as an adjacency list, indexed by ring index.
    graph = {i: set() for i in range(len(carbon_rings))}
    for i in range(len(carbon_rings)):
        for j in range(i+1, len(carbon_rings)):
            if len(carbon_rings[i].intersection(carbon_rings[j])) >= 2:
                graph[i].add(j)
                graph[j].add(i)
    
    # Find connected components of the graph.
    seen = set()
    fused_components = []
    for i in graph:
        if i not in seen:
            # simple DFS
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
            fused_components.append(comp)
    
    # The typical steroid nucleus is 4 fused rings.
    # We look for a component which consists exactly of 4 rings.
    steroid_component = None
    for comp in fused_components:
        if len(comp) == 4:
            steroid_component = comp
            break
    if steroid_component is None:
        return False, ("Fused ring system does not match typical steroid nucleus: "
                       "expected one connected set of 4 carbon-only fused rings, but found fused components of sizes: " +
                       ", ".join(str(len(comp)) for comp in fused_components))
    
    # Gather all atom indices that belong to the steroid nucleus.
    steroid_atom_idxs = set()
    for ring_idx in steroid_component:
        steroid_atom_idxs.update(carbon_rings[ring_idx])
        
    # Check that the fused nucleus has a reasonable number of carbons.
    nucleus_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() in steroid_atom_idxs]
    if len(nucleus_carbons) < 15:
        return False, "Fused ring nucleus has too few carbon atoms to be steroid-like"
    
    # Now look for a suitable ketone group within the steroid nucleus.
    ketone_found = False
    ketone_details = []
    
    # Iterate over atoms in the nucleus that are carbon.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # For our purposes the ketone carbon should be part of at least 2 rings from the steroid nucleus.
        rings_here = [r for r in carbon_rings if atom.GetIdx() in r and r.issubset(steroid_atom_idxs)]
        if len(rings_here) < 2:
            continue  # not situated in two fused rings
        
        for bond in atom.GetBonds():
            # Look for a double bond
            if bond.GetBondTypeAsDouble() != 2:
                continue
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() != 8:
                continue  # must be an oxygen (C=O)
            # Now check that aside from the oxygen, the carbonyl carbon has exactly 2 other neighbors.
            other_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() != other.GetIdx()]
            if len(other_neighbors) != 2:
                continue
            # Both other neighbors must be carbons and must be part of the steroid nucleus.
            if not all(nbr.GetAtomicNum() == 6 and nbr.GetIdx() in steroid_atom_idxs for nbr in other_neighbors):
                continue
            # This appears to be a candidate ketone in the steroid nucleus.
            ketone_found = True
            ketone_details.append(
                f"Carbon {atom.GetIdx()} (in {len(rings_here)} rings) forms a ketone with oxygen {other.GetIdx()}"
            )
            # We only need one correct candidate ketone.
            break
        if ketone_found:
            break

    if not ketone_found:
        return False, "No suitable ring-bound ketone (in the fused steroid nucleus) found – an 11-oxo is required."
    
    reason = ("Molecule contains a steroid-like fused ring system (4 fused carbon-only rings) "
              f"with {len(nucleus_carbons)} nucleus carbons and a candidate ketone group: " +
              "; ".join(ketone_details))
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
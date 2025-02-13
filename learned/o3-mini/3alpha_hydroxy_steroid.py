"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3α-hydroxy steroid
Definition: A 3–hydroxy steroid in which the 3–hydroxy substituent is in the α–position.
This improved version first identifies a fused ring system that is:
  (a) large (at least three rings),
  (b) consisting solely of 5– or 6–membered rings, and
  (c) composed entirely of carbon atoms (as expected for a classical steroid nucleus).
Then it checks that at least one alpha–OH (indicated by the [C@@H](O) SMARTS) is attached
to one of the atoms belonging to that nucleus.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Checks whether the molecule defined by the SMILES string is a 3α-hydroxy steroid.
    Heuristics:
      1. The molecule must contain rings. We get all rings (atom indices) in the molecule.
      2. We build a connectivity graph among rings, connecting two rings if they share two or more atoms.
      3. We select the largest connected ring set (the fused system) and require:
           - At least three rings in the set.
           - Every ring in the set is either 5– or 6–membered.
           - Every atom in the fused system is a carbon (atomic number 6).
      4. We search for a [C@@H](O) substructure which serves as a proxy for an alpha–hydroxyl group,
         and at least one match must reside on an atom that is part of the fused nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a 3α-hydroxy steroid, else False.
        str: Explanation for the decision.
    """
    # Parse SMILES and assign stereochemistry.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found, not a steroid"
    
    # Build a graph of rings: nodes = ring indices; connect two rings if they share >=2 atoms.
    n_rings = len(atom_rings)
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(set(atom_rings[i]).intersection(atom_rings[j])) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # Identify connected components in the ring graph using DFS.
    seen = set()
    fused_components = []
    for i in range(n_rings):
        if i in seen:
            continue
        stack = [i]
        comp = set()
        while stack:
            node = stack.pop()
            if node in comp:
                continue
            comp.add(node)
            for neighbor in adj[node]:
                if neighbor not in comp:
                    stack.append(neighbor)
        seen.update(comp)
        fused_components.append(comp)
    
    if not fused_components:
        return False, "No fused ring systems found, not a steroid nucleus"
    
    # Choose the largest fused ring system (by number of rings).
    largest_component = max(fused_components, key=lambda comp: len(comp))
    if len(largest_component) < 3:
        return False, "Fused ring system too small to be a steroid nucleus"
    
    # Check that every ring in the largest fused set is 5– or 6–membered.
    for i in largest_component:
        ring_size = len(atom_rings[i])
        if ring_size not in (5,6):
            return False, "Fused ring system contains rings not typical for a steroid nucleus (only 5- or 6-membered allowed)"
    
    # Build the set of all atoms that are in the fused ring system.
    fused_atoms = set()
    for idx in largest_component:
        fused_atoms.update(atom_rings[idx])
    
    # Ensure that *every* atom in the fused nucleus is carbon.
    for atom_idx in fused_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return False, "Fused ring nucleus contains heteroatoms; expected all–carbon steroid nucleus"
    
    # Now look for an alpha hydroxyl group using the SMARTS pattern.
    alpha_oh_smarts = "[C@@H](O)"
    alpha_oh_pattern = Chem.MolFromSmarts(alpha_oh_smarts)
    if alpha_oh_pattern is None:
        return False, "Invalid SMARTS for alpha hydroxyl group"
    
    matches = mol.GetSubstructMatches(alpha_oh_pattern)
    if not matches:
        return False, "No [C@@H](O) substructure found (alpha hydroxyl group not detected)"
    
    # Confirm that at least one [C@@H](O) match has its carbon on the fused nucleus.
    for match in matches:
        # match[0] is the carbon bearing the OH.
        if match[0] in fused_atoms:
            return True, "Molecule has a fused (5,6)-ring nucleus composed entirely of carbons and a 3α-hydroxy substituent attached to it"
    
    return False, "Alpha hydroxyl group does not reside on the fused steroid nucleus"

# Testing block:
if __name__ == "__main__":
    # Test with a known example: 5β-Cholestane-3α,7α,12α,23-tetrol
    test_smiles = "C1C[C@H](C[C@]2(C[C@H]([C@@]3([C@](C[C@@H]([C@]4([C@]3(CC[C@@]4([C@@](CC(CC(C)C)O)(C)[H])[H])[H])C)O)([C@@]12C)[H])[H])O)[H])O"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
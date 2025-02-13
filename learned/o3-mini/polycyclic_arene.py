"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: polycyclic arene (polycyclic aromatic hydrocarbon)
A polycyclic aromatic hydrocarbon (PAH) is defined here as a molecule that has 
at least two aromatic rings which are fused (i.e. share one full edge, meaning two atoms)
forming a single connected aromatic substructure. Moreover, the substructure should 
cover a significant fraction (>40%) of all heavy atoms, and be carbon‐dominated (≥60% carbons).
This improved method uses both the overall ring assignments (from GetRingInfo) and a connected
aromatic atom subgraph to capture molecules that otherwise would be missed.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon) 
    based on its SMILES string.
    
    The molecule qualifies if:
      1. It contains at least two fused aromatic rings. 
         Fused here means that rings share a full bond (two atoms) or, if not captured by SSSR,
         a connected aromatic atom subgraph reveals at least two rings.
      2. The fused aromatic substructure covers >40% of the heavy (non-hydrogen) atoms.
      3. At least 60% of the atoms in the fused substructure are carbon.
    
    Args:
      smiles (str): SMILES representation of the molecule.
    
    Returns:
      (bool, str): Tuple where the boolean indicates if the molecule is a PAH,
            and the string gives a reason.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule (this will perceive aromaticity).
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
        
    # Count heavy atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms found in molecule"
    total_heavy = len(heavy_atoms)
    
    # Use the built-in ring information to get aromatic rings that have at least 5 atoms.
    ring_info = mol.GetRingInfo().AtomRings()
    aromatic_rings = []
    for ring in ring_info:
        if len(ring) < 5:
            continue
        # Check that every atom in the ring is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    
    # We first try to detect if at least one connected set of aromatic rings (by shared full edge) exists.
    fused_system = None
    if len(aromatic_rings) >= 2:
        # Build a graph where each node represents one aromatic ring and an edge exists
        # if the rings share at least two atoms (i.e. a full bond).
        ring_graph = {i: set() for i in range(len(aromatic_rings))}
        for i in range(len(aromatic_rings)):
            for j in range(i+1, len(aromatic_rings)):
                common = aromatic_rings[i].intersection(aromatic_rings[j])
                if len(common) >= 2:
                    ring_graph[i].add(j)
                    ring_graph[j].add(i)
        # Find connected ring groups (components) using DFS.
        visited = set()
        fused_components = []
        for i in ring_graph:
            if i not in visited:
                stack = [i]
                comp = set()
                while stack:
                    node = stack.pop()
                    if node in visited:
                        continue
                    visited.add(node)
                    comp.add(node)
                    stack.extend(ring_graph[node] - visited)
                fused_components.append(comp)
        # Take the first component with at least 2 rings.
        for comp in fused_components:
            if len(comp) >= 2:
                # Gather all atom indices in the fused rings.
                fused_atoms = set()
                for idx in comp:
                    fused_atoms.update(aromatic_rings[idx])
                fused_system = fused_atoms
                break

    # If the ring-graph method did not succeed, try an alternative approach:
    # build a graph of aromatic atoms connected by aromatic bonds.
    if fused_system is None:
        aromatic_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
        if not aromatic_indices:
            return False, "No aromatic atoms detected"
        # Build connectivity among aromatic atoms.
        adj = {idx: set() for idx in aromatic_indices}
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetIsAromatic() and a2.GetIsAromatic() and bond.GetIsAromatic():
                idx1 = a1.GetIdx()
                idx2 = a2.GetIdx()
                # Only add if both atoms are in our aromatic list.
                if idx1 in adj and idx2 in adj:
                    adj[idx1].add(idx2)
                    adj[idx2].add(idx1)
        # Find connected components in the aromatic subgraph.
        seen = set()
        components = []
        for idx in aromatic_indices:
            if idx not in seen:
                stack = [idx]
                comp = set()
                while stack:
                    cur = stack.pop()
                    if cur in seen:
                        continue
                    seen.add(cur)
                    comp.add(cur)
                    stack.extend(adj[cur] - seen)
                components.append(comp)
        # Try to select a component that likely corresponds to fused rings.
        # We require that the submolecule from this component contains at least two rings.
        for comp in components:
            # Extract the submolecule from the fused aromatic atoms.
            submol = Chem.PathToSubmol(mol, list(comp))
            # Retrieve its ring info.
            sub_ring_info = submol.GetRingInfo().AtomRings()
            count_arom_rings = 0
            for ring in sub_ring_info:
                if len(ring) >= 5 and all(submol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    count_arom_rings += 1
            if count_arom_rings >= 2:
                fused_system = comp
                break

    if fused_system is None:
        return False, "Aromatic rings are not sufficiently fused (no connected fused system with at least two rings)"
    
    # Now check the fraction of heavy atoms in the fused aromatic system.
    frac_fused = len(fused_system) / float(total_heavy)
    if frac_fused < 0.40:
        return False, f"Fused aromatic system covers only {frac_fused*100:.1f}% of heavy atoms (<40%)"
    
    # Check that the aromatic substructure is largely built from carbons.
    carbons_in_fused = sum(1 for idx in fused_system if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    carbon_frac = carbons_in_fused / float(len(fused_system)) if fused_system else 0
    if carbon_frac < 0.60:
        return False, f"Fused aromatic system is not carbon-dominated ({carbon_frac*100:.1f}% carbons)"
    
    return True, "Molecule is a polycyclic aromatic hydrocarbon with a sufficiently fused aromatic system"


# Example usage for testing (can be removed or commented out if importing this module elsewhere)
if __name__ == "__main__":
    test_smiles = "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34"  # benzo[a]pyrene
    result, reason = is_polycyclic_arene(test_smiles)
    print("SMILES:", test_smiles)
    print("Is polycyclic arene:", result)
    print("Reason:", reason)
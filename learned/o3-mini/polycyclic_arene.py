"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: polycyclic arene (polycyclic aromatic hydrocarbon)
A polycyclic aromatic hydrocarbon (PAH) is defined here as a molecule containing at least 
two aromatic rings that are fused (i.e. share at least one full edge – two or more atoms) 
to form a single connected aromatic substructure. Moreover, the fused system should cover 
a significant fraction (>40%) of the heavy atoms and be largely composed of carbon atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon) based on its SMILES.
    
    The molecule is considered a PAH if:
      1. It contains at least two fully aromatic rings (with 5 or more atoms per ring).
      2. At least two of these rings are sufficiently fused (i.e. they share at least 2 atoms).
      3. The fused aromatic system (the union of atoms in the connected rings) covers 
         a significant fraction (>40%) of all heavy (non-hydrogen) atoms in the molecule.
      4. The fused system is largely composed of carbons (at least 60% of its atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of (True, reason) if molecule qualifies, else (False, reason).
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize (includes aromaticity perception)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # Retrieve all rings in the molecule as sets of atom indices.
    ring_info = mol.GetRingInfo().AtomRings()
    # Filter to rings that are fully aromatic and have at least 5 atoms (to avoid small spurious rings, e.g. 3-membered rings)
    aromatic_rings = []
    for ring in ring_info:
        if len(ring) < 5:
            continue
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    
    if len(aromatic_rings) < 2:
        return False, "Fewer than two fully aromatic rings detected"
    
    # Build a graph where each node is a ring (by its index in aromatic_rings)
    # and add an edge if the rings share at least 2 atoms (i.e., are fused).
    ring_graph = {i: set() for i in range(len(aromatic_rings))}
    for i in range(len(aromatic_rings)):
        for j in range(i+1, len(aromatic_rings)):
            common = aromatic_rings[i].intersection(aromatic_rings[j])
            if len(common) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring graph using DFS.
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
    
    # Look for a fused component that has at least 2 rings.
    fused_system = None
    for comp in fused_components:
        if len(comp) >= 2:
            fused_system = comp
            break
    if fused_system is None:
        return False, "Aromatic rings are not sufficiently fused (no connected fused system with at least two rings)"
    
    # Get union of atom indices that are in the fused rings.
    fused_atoms = set()
    for idx in fused_system:
        fused_atoms.update(aromatic_rings[idx])
    
    # Calculate the fraction of heavy atoms (atomic number > 1) in the entire molecule.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms found"
    frac_fused = len(fused_atoms) / float(len(heavy_atoms))
    
    # We require at least 40% of heavy atoms are in the fused aromatic system.
    if frac_fused < 0.40:
        return False, f"Fused aromatic system covers only {frac_fused*100:.1f}% of heavy atoms (<40%)"
    
    # Optionally, ensure that the fused system is largely composed of carbon.
    carbons_in_fused = sum(1 for idx in fused_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if len(fused_atoms) > 0:
        carbon_frac = carbons_in_fused / float(len(fused_atoms))
    else:
        carbon_frac = 0
    if carbon_frac < 0.60:
        return False, f"Fused aromatic system is not carbon-dominated ({carbon_frac*100:.1f}% carbons)"
    
    return True, "Molecule is a polycyclic aromatic hydrocarbon with a sufficiently fused aromatic system"


# Example usage for testing
if __name__ == "__main__":
    # Example: benzo[a]pyrene, expected to be a PAH.
    test_smiles = "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34"
    result, reason = is_polycyclic_arene(test_smiles)
    print("SMILES:", test_smiles)
    print("Is polycyclic arene:", result)
    print("Reason:", reason)
"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: polycyclic arene (polycyclic aromatic hydrocarbon)

A polycyclic aromatic hydrocarbon (PAH) is defined here as a molecule containing at least 
two aromatic rings that belong to one interconnected aromatic system (i.e. the rings are 
“fused” or connected via aromatic bonds). In addition, the fused aromatic portion of the 
molecule should comprise a significant fraction (here >50%) of the heavy (non‐hydrogen) atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon) from its SMILES.
    
    A molecule is considered a PAH if:
      1. It has at least two aromatic rings (as determined by RDKit routines).
      2. At least two fully aromatic rings appear in one connected aromatic subgraph.
         (Here we build a graph linking aromatic atoms via aromatic bonds.)
      3. The fused aromatic component covers a significant fraction (>50%) of the heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True and a reason if the molecule qualifies, otherwise False and a reason.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Sanitize the molecule (including perception of aromaticity)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # Count the number of aromatic rings in the molecule (using RDKit descriptor)
    n_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_arom_rings < 2:
        return False, "Fewer than two aromatic rings detected"
    
    # Get ring information (atom indices for each ring)
    ring_lists = mol.GetRingInfo().AtomRings()
    # Filter rings to those that are completely aromatic (every atom in the ring is aromatic)
    aromatic_rings = []
    for ring in ring_lists:
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    if len(aromatic_rings) < 2:
        return False, "Fewer than two fully aromatic rings detected"
    
    # Build a graph (as a dict) over aromatic atoms connected by aromatic bonds
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    arom_neighbors = { idx: set() for idx in aromatic_atoms }
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in arom_neighbors and a2 in arom_neighbors:
                arom_neighbors[a1].add(a2)
                arom_neighbors[a2].add(a1)
    
    # Find connected components in the aromatic subgraph using DFS
    visited = set()
    components = []
    for atom in arom_neighbors:
        if atom not in visited:
            stack = [atom]
            comp = set()
            while stack:
                curr = stack.pop()
                if curr in visited:
                    continue
                visited.add(curr)
                comp.add(curr)
                for nbr in arom_neighbors[curr]:
                    if nbr not in visited:
                        stack.append(nbr)
            components.append(comp)
    
    # For each connected component, check how many aromatic rings (from our aromatic_rings list)
    # are completely contained in that component.
    fused_system_found = False
    for comp in components:
        ring_count = 0
        for ring in aromatic_rings:
            if ring.issubset(comp):
                ring_count += 1
        if ring_count >= 2:
            fused_system_found = True
            biggest_arom_comp = comp
            break
            
    if not fused_system_found:
        return False, "Aromatic rings are not sufficiently fused (no connected aromatic system with at least two rings)"
    
    # Check that the fused aromatic system covers a significant fraction of heavy (non-H) atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms found"
    frac_arom = len(biggest_arom_comp) / float(len(heavy_atoms))
    if frac_arom < 0.5:
        return False, f"Fused aromatic system covers only {frac_arom*100:.1f}% of heavy atoms (<50%)"
    
    return True, "Molecule is a polycyclic aromatic hydrocarbon with a fused aromatic system"

# Example usage for testing:
if __name__ == "__main__":
    # Example: benzo[a]pyrene, expected to be True
    test_smiles = "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34"
    result, reason = is_polycyclic_arene(test_smiles)
    print("SMILES:", test_smiles)
    print("Is polycyclic arene:", result)
    print("Reason:", reason)
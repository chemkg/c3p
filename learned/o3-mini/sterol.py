"""
Classifies: CHEBI:15889 sterol
"""
#!/usr/bin/env python
"""
Classifies: Sterol (Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol)
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is defined as any 3-hydroxy steroid whose tetracyclic core is closely related to cholestan-3-ol.
    Its steroid nucleus (gonane) typically consists of 17 carbons in a fused 4-ring system (three six‑membered
    rings and one five‑membered ring) with a hydroxyl group attached to one of the ring carbons.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a sterol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to help detect –OH groups.
    mol = Chem.AddHs(mol)
    
    # Get list of ring atoms (by index)
    ring_atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRing()]
    if not ring_atom_indices:
        return False, "No ring system present, cannot be a sterol"
    
    # Build a connectivity graph for atoms in rings (only bonds connecting ring atoms)
    # Use a dictionary mapping atom index -> set(neighbor indices in ring)
    ring_graph = {idx: set() for idx in ring_atom_indices}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        if idx1 in ring_graph and idx2 in ring_graph:
            ring_graph[idx1].add(idx2)
            ring_graph[idx2].add(idx1)
            
    # Find connected components among ring atoms using DFS.
    visited = set()
    components = []
    for idx in ring_atom_indices:
        if idx not in visited:
            stack = [idx]
            component = set()
            while stack:
                curr = stack.pop()
                if curr not in visited:
                    visited.add(curr)
                    component.add(curr)
                    stack.extend(ring_graph[curr] - visited)
            components.append(component)
    
    # For each ring component, check if it has the steroid nucleus (expect 17 carbon atoms)
    candidate_found = False
    candidate_reason = ""
    for comp in components:
        # Count carbons in this component
        carbon_atoms = [mol.GetAtomWithIdx(i) for i in comp if mol.GetAtomWithIdx(i).GetAtomicNum() == 6]
        if len(carbon_atoms) != 17:
            continue  # Not matching a standard steroid nucleus
        # Now verify that at least one of the ring carbons carries an -OH group.
        hydroxyl_found = False
        for atom in carbon_atoms:
            # Check neighbors outside of the ring component that are oxygen with an attached hydrogen.
            for nbr in atom.GetNeighbors():
                # We expect the hydroxyl oxygen to typically not be in the same ring system.
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in comp:
                    # Check if this oxygen has at least one hydrogen attached.
                    hcount = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                    if hcount >= 1:
                        hydroxyl_found = True
                        break
            if hydroxyl_found:
                break
        if not hydroxyl_found:
            candidate_reason = "Steroid nucleus found but missing a hydroxyl group attached to the ring system"
            continue
        else:
            # Found a candidate steroid nucleus with an -OH
            candidate_found = True
            candidate_reason = "Found a 17-carbon fused ring system with an attached hydroxyl group"
            break

    if candidate_found:
        return True, candidate_reason
    else:
        if candidate_reason:
            return False, candidate_reason
        else:
            return False, "No steroid nucleus (17-carbon fused ring system) identified"

# Example usage (you may remove or comment these out when integrating)
if __name__ == "__main__":
    # Test SMILES for (25R)-3beta,4beta-dihydroxycholest-5-en-26-oic acid (a sterol)
    test_smiles = "C[C@H](CCC[C@@H](C)C(O)=O)[C@H]1CC[C@H]2[C@@H]3CC=C4[C@@H](O)[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_sterol(test_smiles)
    print("Is sterol:", result)
    print("Reason:", reason)
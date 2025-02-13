"""
Classifies: CHEBI:15889 sterol
"""
#!/usr/bin/env python
"""
Classifies: Sterol
Definition: Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol.
The classifier attempts to detect the steroid (gonane) nucleus – a tetracyclic (three six-membered
rings and one five-membered ring, usually containing 17 carbons) and then verifies that at least one of the
core carbons has an attached hydroxyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose tetracyclic core is closely related to cholestan-3-ol.
    This implementation first tries to match a steroid core SMARTS pattern.
    If a match is found, it looks for an -OH (oxygen with at least one hydrogen) attached to one of the core atoms.
    If no SMARTS match is found, it falls back on a ring-analysis approach:
        - It finds fused ring systems (by grouping ring atoms connected via ring bonds).
        - For each connected ring component, it extracts rings (from RDKit ring info) that lie completely
          in the component. It then tries to see if any combination of 4 rings has sizes [5,6,6,6] and the union
          of atoms in the component contains 16-18 carbon atoms.
        - Finally, it checks if one of the atoms in the candidate nucleus carries an external hydroxyl.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a sterol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (helps in detecting -OH groups)
    mol = Chem.AddHs(mol)
    
    # -------- Approach A: SMARTS-based search for a steroid nucleus -----------
    # This SMARTS is chosen to roughly capture the tetracyclic gonane (steroid) core.
    steroid_smarts = "C1CCC2C3CC4CC(C3)C(C4)C2C1"
    steroid_core = Chem.MolFromSmarts(steroid_smarts)
    if steroid_core is None:
        return None, None  # Should not happen
    
    matches = mol.GetSubstructMatches(steroid_core)
    if matches:
        for match in matches:
            # Look for an -OH on any of the core atoms.
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                for nbr in atom.GetNeighbors():
                    # if neighbor is not in the core match and is oxygen
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in match:
                        # Check if oxygen atom is part of an -OH group (has a hydrogen attached)
                        if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                            return True, "Found steroid core (via SMARTS) with an attached hydroxyl group"
        # If we matched the steroid core but none of the atoms carried an external hydroxyl...
        # This captures cases like a deoxy-steroid.
        return False, "Steroid nucleus identified via SMARTS but missing a hydroxyl group on the nucleus"
    
    # -------- Approach B: Fallback using fused ring analysis -----------
    # Get all ring atoms:
    ring_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRing()}
    if not ring_atoms:
        return False, "No ring system present, cannot be a sterol"
    
    # Build a connectivity graph (only including bonds that are in rings) among ring atoms.
    ring_graph = {idx: set() for idx in ring_atoms}
    for bond in mol.GetBonds():
        if bond.IsInRing():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            idx1, idx2 = a1.GetIdx(), a2.GetIdx()
            if idx1 in ring_graph and idx2 in ring_graph:
                ring_graph[idx1].add(idx2)
                ring_graph[idx2].add(idx1)
                
    # Get connected components among ring atoms.
    visited = set()
    fused_components = []
    for idx in ring_atoms:
        if idx not in visited:
            component = set()
            stack = [idx]
            while stack:
                cur = stack.pop()
                if cur not in visited:
                    visited.add(cur)
                    component.add(cur)
                    stack.extend(ring_graph[cur] - visited)
            fused_components.append(component)
    
    # Now use RDKit's ring information:
    all_rings = mol.GetRingInfo().AtomRings()
    
    candidate_found = False
    candidate_reason = ""
    for comp in fused_components:
        # Get atoms in the candidate fused ring system
        comp_atoms = set(comp)
        # Count carbon atoms in the component
        carbons = [mol.GetAtomWithIdx(i) for i in comp_atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 6]
        numC = len(carbons)
        # We expect ~17 carbons in the classical steroid nucleus – allow slight tolerance (16-18)
        if not (16 <= numC <= 18):
            continue
        
        # Get rings completely contained in this component
        comp_rings = [r for r in all_rings if set(r).issubset(comp_atoms)]
        if len(comp_rings) < 4:
            continue  # need 4 rings in a steroid nucleus
        
        # Check if any combination of 4 rings has sizes that when sorted equals [5,6,6,6]
        valid_combo = False
        for ring_combo in itertools.combinations(comp_rings, 4):
            sizes = sorted([len(ring) for ring in ring_combo])
            if sizes == [5, 6, 6, 6]:
                valid_combo = True
                break
        if not valid_combo:
            continue
        
        # Look for hydroxyl group attached (oxygen with hydrogen) to any atom in this component,
        # but the oxygen should lie outside of the fused ring system.
        hydroxyl_found = False
        for idx in comp_atoms:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbons in the nucleus
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in comp_atoms:
                    # if at least one hydrogen is attached to the oxygen, treat it as an -OH
                    if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                        hydroxyl_found = True
                        break
            if hydroxyl_found:
                break
        
        if not hydroxyl_found:
            candidate_reason = "Steroid nucleus found in fused ring system but missing a hydroxyl group attached to the nucleus"
            continue
        else:
            candidate_found = True
            candidate_reason = "Found fused steroid nucleus (4 rings with sizes 5,6,6,6 and ~17 carbons) with an attached hydroxyl group"
            break
    
    if candidate_found:
        return True, candidate_reason
    else:
        if candidate_reason:
            return False, candidate_reason
        else:
            return False, "No steroid nucleus identified"

# Example usage:
if __name__ == "__main__":
    # Test with one provided sterol SMILES: (25R)-3beta,4beta-dihydroxycholest-5-en-26-oic acid
    test_smiles = "C[C@H](CCC[C@@H](C)C(O)=O)[C@H]1CC[C@H]2[C@@H]3CC=C4[C@@H](O)[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    res, reason = is_sterol(test_smiles)
    print("Is sterol:", res)
    print("Reason:", reason)
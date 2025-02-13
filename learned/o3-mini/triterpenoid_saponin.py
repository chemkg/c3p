"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
#!/usr/bin/env python
"""
Classifies: Triterpenoid saponin
A triterpenoid saponin is defined as a terpene glycoside in which the terpene (aglycone)
moiety is a triterpenoid. That is, the molecule must contain one or more sugar units (detected as sugar rings)
attached via O–glycosidic bonds to an aglycone that has roughly 30 carbon atoms (with some margin)
and a fused ring system (typically four or more rings) and a high sp3 carbon fraction.
This heuristic approach uses substructure detection of candidate sugar rings and counts the aglycone atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdchem

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    
    Checks performed:
      1. Parses the SMILES.
      2. Looks for at least one sugar (glycoside) ring. Here, a candidate sugar ring is
         defined heuristically as a 5- or 6-membered ring that contains exactly one oxygen
         atom (typical for pyranose or furanose rings) and that has at least two hydroxyl (–OH)
         substituents on ring carbon(s).
      3. Defines the aglycone as the set of all atoms not assigned to a detected sugar ring.
      4. Counts the number of carbon atoms in the aglycone and requires that this number be
         roughly around 30 (we allow from 25 to 40).
      5. Counts the rings that are completely contained in the aglycone; at least 4 rings are expected.
      6. As an extra check for a terpene (triterpenoid) scaffold, we compute the fraction of sp3 carbons
         among aglycone carbons; triterpenoids are largely saturated so we require a minimum sp3 fraction.
         
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a triterpenoid saponin, False otherwise.
        str: A reason describing the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    
    # Heuristic to detect sugar rings: 5- or 6-membered ring with exactly 1 O atom and at least 2 hydroxyl substituents.
    sugar_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        # Count oxygen atoms in ring
        num_oxygen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if num_oxygen != 1:
            continue
        
        # Count hydroxyl substituents attached to ring carbons.
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # carbon atom
                # A hydroxyl is an O neighbor having degree 1.
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                        hydroxyl_count += 1
        if hydroxyl_count >= 2:
            sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar moiety (glycoside unit) detected"
    
    # Union all detected sugar ring atoms.
    sugar_atoms = set()
    for ring_set in sugar_rings:
        sugar_atoms |= ring_set
    
    # Define aglycone as atoms not in any sugar ring.
    total_atoms = set(range(mol.GetNumAtoms()))
    aglycone_atoms = total_atoms - sugar_atoms
    
    # Count carbon atoms in the aglycone.
    aglycone_carbons = sum(1 for idx in aglycone_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Triterpenoid aglycones are typically built from ~30 carbons; allow a range (e.g. 25 to 40).
    if not (25 <= aglycone_carbons <= 40):
        return False, f"Aglycone carbon count ({aglycone_carbons}) not within expected range for triterpenoid (~30 carbons)"
    
    # Count rings entirely within the aglycone.
    aglycone_ring_count = 0
    for ring in ring_info.AtomRings():
        if set(ring).issubset(aglycone_atoms):
            aglycone_ring_count += 1
    if aglycone_ring_count < 4:
        return False, f"Fused ring system in aglycone insufficient (found {aglycone_ring_count} rings, expected at least 4)"
    
    # As an extra check for terpene-like saturation, compute sp3 carbon ratio.
    aglycone_sp3 = 0
    for idx in aglycone_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == rdchem.HybridizationType.SP3:
            aglycone_sp3 += 1
    sp3_ratio = aglycone_sp3 / aglycone_carbons if aglycone_carbons > 0 else 0
    if sp3_ratio < 0.6:
        return False, f"Aglycone sp3 carbon fraction too low ({sp3_ratio:.2f}, expected > 0.6)"
    
    return True, "Molecule contains a triterpenoid aglycone with appropriate carbon count and fused rings plus a sugar moiety"

# Example usage; you may remove or adjust the example below.
if __name__ == "__main__":
    # Example SMILES from the prompt (ginsenoside Re)
    test_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_triterpenoid_saponin(test_smiles)
    print(result, ":", reason)
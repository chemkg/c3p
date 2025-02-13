"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
#!/usr/bin/env python
"""
Classifies: Triterpenoid saponin
A triterpenoid saponin is defined as a terpene glycoside in which the terpene (aglycone)
moiety is a triterpenoid. That is, the molecule must contain one or more sugar units attached
via an O–glycosidic bond to an aglycone that is largely based on a ~30-carbon fused ring system.
This heuristic approach first searches for candidate sugar rings (5- or 6-membered with one ring oxygen
and an O-connection to the rest of the molecule), then defines the aglycone as atoms not included in
any sugar unit. It then checks that the aglycone has a carbon count in a broadened range, enough fused rings,
and a sufficient sp3 carbon fraction.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdchem

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES.
    
    Heuristic checks performed:
      1. Parse the SMILES.
      2. Detect sugar rings:
           - Consider only rings of size 5 or 6.
           - Require exactly one oxygen atom in the ring.
           - The ring must be attached (via an oxygen bond) to an atom outside the ring.
      3. Define the aglycone as the atoms not contained in any detected sugar ring.
      4. Count carbon atoms in the aglycone; expect roughly 25–60 carbons.
      5. Count rings that are completely contained in the aglycone; require at least 4.
      6. Compute the sp3 fraction among aglycone carbons (require at least 50%).
      
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a triterpenoid saponin, False otherwise.
        str: A reason message describing the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    
    # Step 2: Detect candidate sugar rings.
    sugar_rings = []
    # Loop over all rings in the molecule.
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        # Count the number of oxygen atoms in the ring.
        num_ring_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if num_ring_oxygens != 1:
            continue
        # Check if at least one atom in the ring has an oxygen neighbor outside the ring.
        attached_via_oxygen = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # If neighbor is not in the ring and is oxygen, mark the ring as glycosidically attached.
                if nbr_idx not in ring and nbr.GetAtomicNum() == 8:
                    attached_via_oxygen = True
                    break
            if attached_via_oxygen:
                break
        if attached_via_oxygen:
            sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar moiety (glycoside unit) detected"
    
    # Step 3: Define aglycone as all atoms not found in any sugar ring.
    sugar_atoms = set()
    for ring_set in sugar_rings:
        sugar_atoms |= ring_set
    total_atoms = set(range(mol.GetNumAtoms()))
    aglycone_atoms = total_atoms - sugar_atoms
    
    # Step 4: Count carbon atoms in the aglycone.
    aglycone_carbons = sum(1 for idx in aglycone_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if not (25 <= aglycone_carbons <= 60):
        return False, f"Aglycone carbon count ({aglycone_carbons}) not within expected range for triterpenoid (~25–60 carbons)"
    
    # Step 5: Count rings entirely contained within the aglycone.
    aglycone_ring_count = 0
    for ring in ring_info.AtomRings():
        if set(ring).issubset(aglycone_atoms):
            aglycone_ring_count += 1
    if aglycone_ring_count < 4:
        return False, f"Fused ring system in aglycone insufficient (found {aglycone_ring_count} rings, expected at least 4)"
    
    # Step 6: Check aglycone saturation by computing sp3 fraction among aglycone carbons.
    aglycone_sp3 = 0
    for idx in aglycone_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == rdchem.HybridizationType.SP3:
            aglycone_sp3 += 1
    sp3_ratio = aglycone_sp3 / aglycone_carbons if aglycone_carbons > 0 else 0
    if sp3_ratio < 0.50:
        return False, f"Aglycone sp3 carbon fraction too low ({sp3_ratio:.2f}, expected > 0.50)"
    
    return True, "Molecule contains a triterpenoid aglycone with appropriate carbon count, fused ring system, and a sugar moiety"

# Example usage; you may adjust or remove depending on your workflow.
if __name__ == "__main__":
    # Example: ginsenoside Re (a known triterpenoid saponin)
    test_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_triterpenoid_saponin(test_smiles)
    print(result, ":", reason)
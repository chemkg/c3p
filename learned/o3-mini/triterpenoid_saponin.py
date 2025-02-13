"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
#!/usr/bin/env python
"""
Classifies: Triterpenoid saponin

A triterpenoid saponin is defined as a terpene glycoside in which the terpene (aglycone)
moiety is a triterpenoid. This heuristic approach:
  1. Parses the SMILES into an RDKit molecule.
  2. Detects candidate sugar rings:
       - Considers rings of size 5 or 6.
       - Requires exactly one oxygen atom (the ring heteroatom) in the ring.
       - Requires that at the ring’s anomeric (or any) carbon at least one exocyclic oxygen
         is bonded to a carbon outside the ring.
  3. Defines the aglycone as the atoms not assigned to any sugar ring.
  4. Checks that the aglycone contains a reasonable number of carbon atoms (roughly 25–60).
  5. Groups aglycone rings into fused clusters (rings sharing ≥2 atoms) and requires
     that at least one cluster contains three or more rings.
  6. Computes the fraction of sp3 hybridized carbons within the aglycone (target ≥0.55).
  
Because the nature of triterpenoid saponins is subtle, this is but one heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdchem

def get_fused_ring_clusters(ring_list):
    """
    Given a list of rings (each a set of atom indices) create a list of fused ring clusters.
    Two rings are considered fused if they share ≥2 atoms.
    We form clusters by merging rings that are connected.
    """
    clusters = []
    for ring in ring_list:
        # Try to add this ring to an existing cluster if it shares enough atoms.
        added = False
        for cluster in clusters:
            # If ring fuses with any ring in the cluster, merge.
            if any(len(ring & other) >= 2 for other in cluster):
                cluster.append(ring)
                added = True
                break
        if not added:
            clusters.append([ring])
    return clusters

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES.
    
    Heuristic checks performed:
      1. Parse the SMILES.
      2. Detect sugar rings:
           - Only consider rings of size 5 or 6.
           - The ring must include exactly one oxygen atom.
           - At least one atom in the ring should be connected via an oxygen to an atom
             outside the ring (and that oxygen must be connected to at least one carbon).
      3. Define aglycone as atoms not in any sugar ring.
      4. Check that the aglycone has a reasonable carbon count (25–60).
      5. Identify fused ring clusters in the aglycone.
         Require that at least one cluster has ≥3 rings.
      6. Compute the sp3 fraction among aglycone carbons (require at least 55%).
      
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a triterpenoid saponin,
              False otherwise.
        str: A reason message describing the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    
    # Step 2: Identify candidate sugar rings.
    sugar_rings = []
    for ring in all_rings:
        if len(ring) not in [5,6]:
            continue
        # Count oxygen atoms in this ring.
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        oxy_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        if oxy_count != 1:
            continue
        # Check for an exocyclic oxygen bond.
        attached_through_oxygen = False
        for i in ring:
            atom = mol.GetAtomWithIdx(i)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring:
                    continue
                # Check if the connecting bond is via oxygen:
                # Either the neighbor itself is O and it is bonded to at least one carbon
                if nbr.GetAtomicNum() == 8:
                    # make sure that neighbor is attached to a non-sugar carbon (tentative)
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() not in ring and nn.GetAtomicNum() == 6:
                            attached_through_oxygen = True
                            break
                if attached_through_oxygen:
                    break
            if attached_through_oxygen:
                break
        if attached_through_oxygen:
            sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar moiety (glycoside unit) detected"
    
    # Step 3: Define aglycone as atoms not in any detected sugar ring.
    sugar_atom_idxs = set()
    for ring_set in sugar_rings:
        sugar_atom_idxs |= ring_set
    total_idxs = set(range(mol.GetNumAtoms()))
    aglycone_idxs = total_idxs - sugar_atom_idxs
    if not aglycone_idxs:
        return False, "No aglycone atoms remaining after removing sugar rings"

    # Step 4: Count carbon atoms in aglycone.
    aglycone_carbons = sum(1 for idx in aglycone_idxs if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if not (25 <= aglycone_carbons <= 60):
        return False, f"Aglycone carbon count ({aglycone_carbons}) not within expected range for triterpenoid (~25–60 carbons)"
    
    # Step 5: Count fused rings within the aglycone.
    # Collect rings that lie entirely in the aglycone.
    aglycone_rings = []
    for ring in all_rings:
        if set(ring).issubset(aglycone_idxs):
            aglycone_rings.append(set(ring))
    if not aglycone_rings:
        return False, "No rings found in aglycone"
    
    clusters = get_fused_ring_clusters(aglycone_rings)
    largest_cluster_size = max(len(cluster) for cluster in clusters) if clusters else 0
    # We require at least 3 fused rings in one cluster.
    if largest_cluster_size < 3:
        return False, f"Fused ring system in aglycone insufficient (largest cluster has {largest_cluster_size} rings, expected ≥3)"
    
    # Step 6: Compute the sp3 fraction among aglycone carbons.
    sp3_count = 0
    for idx in aglycone_idxs:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == rdchem.HybridizationType.SP3:
            sp3_count += 1
    sp3_fraction = sp3_count / aglycone_carbons if aglycone_carbons > 0 else 0
    if sp3_fraction < 0.55:
        return False, f"Aglycone sp3 carbon fraction too low ({sp3_fraction:.2f}, expected ≥0.55)"
    
    return True, "Molecule contains a triterpenoid aglycone with appropriate carbon count, fused ring system, and a sugar moiety"

# Example usage:
if __name__ == "__main__":
    # Test with one known triterpenoid saponin (ginsenoside Re) SMILES.
    test_smiles = "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_triterpenoid_saponin(test_smiles)
    print(result, ":", reason)
"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Naturally occurring prostaglandins (derived from C20 prostanoic acid)

A prostaglandin is heuristically defined here as a molecule that:
  1. Contains at least 20 carbon atoms.
  2. Contains at least one non‐aromatic cyclopentane (5‐membered) ring.
  3. Contains a carboxyl (or ester) group – defined by a C(=O)OX fragment – that is “anchored” 
     to the cyclopentane ring (i.e. the carboxyl carbon is at most 2 bonds from one of the ring atoms).
  4. Contains a long aliphatic chain – here at least 4 connected sp3 carbons – that is attached 
     directly (by at least one bond) to the cyclopentane ring.
  5. The cyclopentane ring must have at least two distinct attachment sites (one for the carboxyl/ester 
     and one for the chain).
    
Note: This heuristic is an approximation. Some prostaglandin derivatives may have extra linkers 
or differing substitution. Also, some non-prostaglandin molecules may meet these criteria.
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule qualifies as a prostaglandin (or prostaglandin derivative) based
    on a heuristic. It uses the following criteria:
      1. At least 20 carbon atoms.
      2. At least one non-aromatic cyclopentane ring.
      3. A carboxyl/ester group (C(=O)OX) whose carbon is within 2 bonds (graph distance) 
         of a cyclopentane ring atom.
      4. At least one long aliphatic chain (4 connected sp3 carbons) that is directly 
         attached to the cyclopentane ring.
      5. The cyclopentane ring displays at least two distinct attachments (one from the carboxyl
         group and one from the chain) coming off different ring atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the prostaglandin criteria, False otherwise.
        str: Reason message for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for at least 20 carbons.
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 20:
        return False, f"Insufficient carbons: found {n_carbons} but need at least 20 for a prostaglandin"
    
    # 2. Identify non-aromatic cyclopentane rings.
    ring_info = mol.GetRingInfo().AtomRings()
    cyclopentane_rings = []
    for ring in ring_info:
        if len(ring) == 5 and all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            cyclopentane_rings.append(set(ring))
    if not cyclopentane_rings:
        return False, "No non-aromatic cyclopentane (5-membered) ring found, which is required for a prostaglandin scaffold"
    
    # 3. Identify carboxyl (or ester) groups.
    # This pattern matches a C(=O)O fragment where the oxygen can be -OH, -OR, or even anionic.
    carboxyl_smarts = "[CX3](=O)[OX2H1,O-]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl/ester group (C(=O)OX) found, which is required for prostanoic acid derivatives"
    
    # 4. Identify a long aliphatic chain composed of at least 4 connected sp3 (tetrahedral) carbons.
    chain_smarts = "[CX4]-[CX4]-[CX4]-[CX4]"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No aliphatic chain (at least 4 connected sp3 carbons) detected"
    
    # Precompute distance matrix to measure bond distances.
    dmat = Chem.GetDistanceMatrix(mol)
    
    # For each cyclopentane ring, we will see if two separate attachment sites exist:
    # one for a nearby carboxyl group and one for an attached aliphatic chain.
    qualifying_ring_found = False
    for ring in cyclopentane_rings:
        ring = set(ring)  # for easier membership testing
        # These sets will collect ring atoms (by index) that are connected (within the desired distance)
        # to a carboxyl carbon and to a chain fragment.
        carboxyl_attached_atoms = set()
        chain_attached_atoms = set()
        
        # Check carboxyl: for each carboxyl match, take the carbon (first atom in the match) 
        # and see if it lies within bond distance 2 of any ring atom.
        for match in carboxyl_matches:
            carboxyl_c = match[0]  # the carbonyl carbon index
            for aidx in ring:
                # distance 1 means directly attached; distance 2 means one bridging atom.
                if dmat[aidx][carboxyl_c] <= 2:
                    carboxyl_attached_atoms.add(aidx)
        
        # Check the chain attachment: for each chain match, see if any atom of the chain is a neighbor 
        # of a ring atom.
        for chain_match in chain_matches:
            chain_atoms = set(chain_match)
            for aidx in ring:
                atom = mol.GetAtomWithIdx(aidx)
                # Look at neighbors of this ring atom that are not in the ring.
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in ring and nbr_idx in chain_atoms:
                        chain_attached_atoms.add(aidx)
                        break  # one attachment is enough for this ring atom
        # For a valid prostaglandin scaffold, the cyclopentane ring should have at least 2 distinct attachment sites;
        # one for the carboxyl group and one for the chain. They need not be on different atoms if a single atom has both 
        # but here we require that the union of attachment sites has size >= 2.
        if carboxyl_attached_atoms and chain_attached_atoms:
            if len(carboxyl_attached_atoms.union(chain_attached_atoms)) >= 2:
                qualifying_ring_found = True
                break
    
    if not qualifying_ring_found:
        # Decide on a reason based on whether the carboxyl or chain wasn’t attached closely
        if not carboxyl_matches:
            return False, "No carboxyl group found to anchor the prostaglandin scaffold"
        else:
            return False, "Found a carboxyl/ester group and a long aliphatic chain, but they are not appropriately anchored on the cyclopentane ring"
    
    return True, ("Molecule meets prostaglandin criteria: at least 20 carbons, a non-aromatic cyclopentane ring with a carboxyl/ester "
                  "group (directly or via one bridging bond) and an attached long aliphatic chain (with at least 4 connected sp3 carbons) "
                  "coming off distinct ring atoms.")

# Example usage:
# Uncomment the lines below to test one of the prostaglandin SMILES.
# Example: prostaglandin E1
# smiles_example = "CCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1CCCCCCC(O)=O"
# result, reason = is_prostaglandin(smiles_example)
# print(result, reason)
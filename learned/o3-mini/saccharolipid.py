"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI saccharolipid (Lipids that contain a carbohydrate moiety.)
This version uses two strategies:
  1. Identify candidate sugar rings: non‐aromatic 5– or 6–membered rings made solely of C and O (with at least one O).
  2. For each sugar ring, examine each non‐ring neighbor. Either: 
      a. if it is an oxygen that forms an ester bond (bonded to a carbonyl carbon), 
         then follow the acyl (fatty) chain off that carbon to check for a long (>=7 C) alkyl chain, or 
      b. if directly attached to a carbon (outside the ring) that is part of an aliphatic chain, follow that chain.
  Only if one such connection is found will the molecule be classified as a saccharolipid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid that contains a carbohydrate moiety.
    Here we require that the molecule contains at least one sugar ring—specifically,
    a non‐aromatic 5– or 6–membered ring containing only C and O (and at least one O)—
    and that at least one of these sugar rings is covalently linked (directly or via an ester/amide bond)
    to a long, aliphatic chain (heuristically, one showing at least 7 contiguous aliphatic carbon atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if it is a saccharolipid, False otherwise.
        str: A reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # -------------------------------
    # Step 1. Identify candidate sugar rings.
    # Look for rings with size 5 or 6, check that every atom in the ring is either C (6) or O (8),
    # require at least one O, and that none of these atoms is aromatic.
    ring_info = mol.GetRingInfo()
    sugar_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        atomic_nums = [atom.GetAtomicNum() for atom in atoms]
        # must contain only C or O and at least one O
        if all(num in (6, 8) for num in atomic_nums) and any(num == 8 for num in atomic_nums):
            # also require that none of these atoms is aromatic
            if not any(atom.GetIsAromatic() for atom in atoms):
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety (5- or 6-membered, non-aromatic ring of C/O) detected"
    
    # -------------------------------
    # Helper: Given a starting atom (typically the carbonyl carbon) try to follow an acyclic aliphatic chain.
    # We use a depth-first search over neighbors that are sp3 carbons (atomic num 6) not in any ring,
    # and count how many such carbons in a contiguous chain.
    def follow_aliphatic_chain(start_idx, visited=None):
        if visited is None:
            visited = set()
        count = 0
        stack = [start_idx]
        while stack:
            idx = stack.pop()
            if idx in visited:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            # check: should be sp3 carbon, not aromatic, and not in a ring
            if atom.GetAtomicNum() == 6 and (atom.GetHybridization().name == "SP3") \
                 and (not atom.IsInRing()) and (not atom.GetIsAromatic()):
                count += 1
                for nb in atom.GetNeighbors():
                    nb_idx = nb.GetIdx()
                    # Only follow bonds between carbons
                    if nb.GetAtomicNum() == 6 and nb_idx not in visited:
                        stack.append(nb_idx)
        return count

    # -------------------------------
    # Step 2. For each sugar ring, look at bonds from ring atoms to outside atoms.
    # Try to detect attachment to a long acyl chain.
    found_attachment = False
    for ring in sugar_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                # Ignore if neighbor is also in the sugar ring.
                if nb_idx in ring:
                    continue
                # -- Strategy A: Look if the sugar is attached via an oxygen to an ester-like linkage.
                if nb.GetAtomicNum() == 8 and (not nb.GetIsAromatic()):
                    # For the neighbor oxygen, check its other neighbors.
                    for oxygen_nb in nb.GetNeighbors():
                        # Skip going back to the sugar ring atom.
                        if oxygen_nb.GetIdx() == idx:
                            continue
                        # Look for a carbon with a double bond to oxygen (a carbonyl).
                        if oxygen_nb.GetAtomicNum() == 6:
                            has_carbonyl = False
                            for candidate in oxygen_nb.GetNeighbors():
                                # check bond order; we use GetExplicitValence() is not precise so instead check symbol
                                # Here a carbonyl oxygen usually is represented as =O, so we look for a neighbor with atomic num 8
                                if candidate.GetAtomicNum() == 8:
                                    # We test bond type:
                                    bond = mol.GetBondBetweenAtoms(oxygen_nb.GetIdx(), candidate.GetIdx())
                                    if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                                        has_carbonyl = True
                                        break
                            if has_carbonyl:
                                # Now, follow the chain starting from the carbonyl carbon (oxygen_nb)
                                chain_length = follow_aliphatic_chain(oxygen_nb.GetIdx())
                                if chain_length >= 7:
                                    found_attachment = True
                                    break
                    if found_attachment:
                        break
                # -- Strategy B: Direct C–C link: if the sugar ring carbon is bonded to an external aliphatic carbon.
                if nb.GetAtomicNum() == 6 and (not nb.GetIsAromatic()) and (not nb.IsInRing()):
                    # follow the chain starting at nb.
                    chain_length = follow_aliphatic_chain(nb.GetIdx())
                    if chain_length >= 7:
                        found_attachment = True
                        break
            if found_attachment:
                break
        if found_attachment:
            break

    if not found_attachment:
        return False, "Sugar (carbohydrate) moiety not covalently attached to a long acyl chain (>=7 contiguous aliphatic carbons) by ester or direct C–C bond"

    # Optionally, one could add further criteria (e.g. relative molecular weight) but we stop here.
    return True, "Molecule contains a carbohydrate ring attached to a long acyl (lipid) chain"

# Example usage
if __name__ == "__main__":
    # Example saccharolipid (one of the provided examples)
    smiles_example = "CCCCCCCCCCCCC(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC"
    result, reason = is_saccharolipid(smiles_example)
    print("Is saccharolipid?", result)
    print("Reason:", reason)
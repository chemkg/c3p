"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI saccharolipid (Lipids that contain a carbohydrate moiety.)

This version improves on the previous strategy by:
  1. Defining a sugar ring as a non‐aromatic 5– or 6–membered ring composed solely of C and O (and containing at least one O).
  2. For each sugar ring, examining all bonds (from ring atoms to external atoms) to see if the ring is covalently attached to a long alkyl (acyl) chain.
     – Strategy A: if attached via an oxygen (as in an ester bond sugar–O–C(=O)–R) then follow from the carbonyl carbon’s non‐carbonyl neighbor.
     – Strategy B: if attached directly (a C–C bond) then follow that chain.
     
A helper function “longest_alkyl_chain” is defined that obeys strict criteria:
   • Only sp3 (non‐aromatic) carbon atoms not belonging to any ring are included.
   • The search is limited to a contiguous segment (ignoring branches) using recursion.

If any sugar ring is found to have an attachment (via either strategy) to an acyl chain of at least 7 contiguous carbon atoms, the molecule is classified as a saccharolipid.
"""

from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid with a carbohydrate (sugar) moiety that is covalently attached
    (by an ester or direct C–C bond) to a long aliphatic chain (>= 7 contiguous sp3 non‐aromatic carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if saccharolipid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------------------------------------------------------------------
    # Step 1. Identify candidate sugar rings.
    # Our sugar ring is defined as a 5- or 6-membered ring, all atoms are C or O,
    # at least one oxygen is present, and no atom is aromatic.
    ring_info = mol.GetRingInfo()
    sugar_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Accept only rings with atoms that are carbon (6) or oxygen (8),
        # and require at least one oxygen.
        atomic_nums = [atom.GetAtomicNum() for atom in atoms]
        if all(num in (6,8) for num in atomic_nums) and any(num == 8 for num in atomic_nums):
            if not any(atom.GetIsAromatic() for atom in atoms):
                # Save as a set for fast membership test
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety (non‐aromatic 5- or 6-membered ring of C/O) detected"
    
    # ---------------------------------------------------------------------
    # Helper: recursively find the longest chain of contiguous sp3, non-aromatic, non-ring carbons.
    def longest_alkyl_chain(atom_idx, exclude, visited=None):
        """
        Returns the length (number of carbons) of the longest contiguous chain starting
        from the given atom (which should be an sp3 carbon, non-ring) and not traversing atoms in 'exclude'.
        """
        if visited is None:
            visited = set()
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        max_length = 1  # count the starting atom
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            # Only follow if neighbor is carbon, not aromatic, sp3 and not in any ring.
            if nb_idx in visited or nb_idx in exclude:
                continue
            if (nb.GetAtomicNum() == 6 and 
                not nb.GetIsAromatic() and 
                nb.GetHybridization().name == "SP3" and
                not nb.IsInRing()):
                length = 1 + longest_alkyl_chain(nb_idx, exclude, visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    MIN_CHAIN_LENGTH = 7  # the threshold for a long acyl chain

    # ---------------------------------------------------------------------
    # Strategy evaluation: For each sugar ring, check all bonds from a ring atom to an external atom.
    found_attachment = False
    for ring in sugar_rings:
        for idx in ring:
            sugar_atom = mol.GetAtomWithIdx(idx)
            for nb in sugar_atom.GetNeighbors():
                nb_idx = nb.GetIdx()
                if nb_idx in ring:
                    continue  # neighbor is part of the sugar ring
                
                # Strategy A: Attachment via oxygen (as in an ester bond sugar–O–C(=O)–R)
                if nb.GetAtomicNum() == 8 and (not nb.GetIsAromatic()):
                    # nb is an oxygen bridging the sugar ring.
                    # Look at its neighbors (besides the sugar_atom) to see if one is a carbonyl carbon.
                    for o_nb in nb.GetNeighbors():
                        if o_nb.GetIdx() == idx:
                            continue
                        if o_nb.GetAtomicNum() == 6:
                            # Check if o_nb (the candidate carbonyl carbon) has a double-bonded oxygen.
                            has_carbonyl = False
                            for c_nb in o_nb.GetNeighbors():
                                # Look for an oxygen with a double bond to o_nb.
                                if c_nb.GetAtomicNum() == 8:
                                    bond = mol.GetBondBetweenAtoms(o_nb.GetIdx(), c_nb.GetIdx())
                                    if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                                        has_carbonyl = True
                                        break
                            if has_carbonyl:
                                # Now (in a typical ester) the acyl chain comes off o_nb (but not via its double-bonded oxygens).
                                # So check each neighbor of o_nb other than nb and sugar_atom:
                                for acyl_candidate in o_nb.GetNeighbors():
                                    candidate_idx = acyl_candidate.GetIdx()
                                    if candidate_idx in (nb.GetIdx(), idx):
                                        continue
                                    # We require candidate to be a carbon with sp3 geometry (i.e. part of the alkyl chain).
                                    if (acyl_candidate.GetAtomicNum() == 6 and 
                                        not acyl_candidate.GetIsAromatic() and 
                                        acyl_candidate.GetHybridization().name == "SP3" and
                                        not acyl_candidate.IsInRing()):
                                        chain_length = longest_alkyl_chain(candidate_idx, exclude=ring)
                                        if chain_length >= MIN_CHAIN_LENGTH:
                                            found_attachment = True
                                            break
                                if found_attachment:
                                    break
                    if found_attachment:
                        break
                # Strategy B: Direct C–C bond from sugar to an external aliphatic carbon.
                if nb.GetAtomicNum() == 6 and (not nb.GetIsAromatic()) and (nb.GetHybridization().name == "SP3") and (not nb.IsInRing()):
                    chain_length = longest_alkyl_chain(nb_idx, exclude=ring)
                    if chain_length >= MIN_CHAIN_LENGTH:
                        found_attachment = True
                        break
            if found_attachment:
                break
        if found_attachment:
            break

    if not found_attachment:
        return False, "Sugar moiety not covalently attached to a long acyl chain (>=7 contiguous sp3 carbons) by ester or direct C–C bond"
    else:
        return True, "Molecule contains a carbohydrate ring attached to a long acyl (lipid) chain"

# Example usage:
if __name__ == "__main__":
    # One example saccharolipid (taken from supplied examples)
    example_smiles = "CCCCCCCCCCCCC(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(O)(=O)=O)[C@@H]1OC(=O)CCCCCCCCCCCCCCC"
    result, reason = is_saccharolipid(example_smiles)
    print("Is saccharolipid?", result)
    print("Reason:", reason)
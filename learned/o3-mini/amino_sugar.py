"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies amino sugars:
Definition: An amino sugar is a sugar molecule in which one or more alcoholic hydroxyl groups attached to the sugar‐ring 
             have been replaced by unsubstituted or substituted amino groups.
             
The improved strategy:
  1. Parse the SMILES string.
  2. Identify candidate sugar rings. We require:
       • ring size of 5 or 6 atoms 
       • non-aromatic ring 
       • exactly one oxygen atom in the ring (typical for pyranoses or furanoses)
  3. For each candidate ring, inspect the exocyclic substituents from ring carbons:
       • Count –OH groups (an oxygen not in the ring bonded to at least one hydrogen)
       • Look for at least one free –NH2 group (a nitrogen not in the ring that carries two explicit hydrogens and is not part of an amide)
       • To qualify as a sugar (for our purposes) the candidate ring should have several (≥3) polar substituents.
  4. Return True if at least one candidate ring meets these criteria.
  
If no suitable candidate ring is found or no valid amino substituent is found on a sugar ring, return False.
"""

from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an amino sugar, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information from the molecule
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    
    # Identify candidate rings: (5- or 6-membered, only one oxygen in the ring, non-aromatic)
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        # Retrieve ring atoms and check aromaticity.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue  # skip aromatic rings
        # Count oxygen atoms in this ring
        oxygen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        if oxygen_count != 1:
            continue
        # Save candidate ring (as a set for fast membership testing)
        candidate_rings.append(set(ring))
    
    if not candidate_rings:
        return False, "No sugar-like ring (5- or 6-membered, non-aromatic with exactly one ring oxygen) detected"
    
    # Define helper function to check if an exocyclic substituent is an -OH group.
    def is_hydroxyl(neighbor):
        # Check if the substituent atom is oxygen, has at least one hydrogen, and is not in any ring.
        if neighbor.GetAtomicNum() != 8:
            return False
        # Count explicit hydrogens. For our purpose, at least one explicit hydrogen.
        numH = neighbor.GetTotalNumHs(includeNeighbors=True)
        return numH >= 1  # assume OH if bonded to H (we could further restrict the degree)

    # Define helper function to check if an exocyclic substituent is a free -NH2 group.
    def is_free_amino(neighbor):
        if neighbor.GetAtomicNum() != 7:
            return False
        # Check that the nitrogen has exactly two hydrogens (or two total H) as substituents.
        numH = neighbor.GetTotalNumHs(includeNeighbors=True)
        if numH != 2:
            return False
        # Also, check that none of its bonds is a double bond to an oxygen (which would indicate an amide)
        for bond in neighbor.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(neighbor)
                if other.GetAtomicNum() == 8:
                    return False
        return True

    # For each candidate ring, check exocyclic substituents on ring carbons.
    for ring_set in candidate_rings:
        polar_substituent_count = 0
        found_amino = False
        # Iterate over atoms in the candidate ring. Sugar rings are mostly made up of carbons except the one oxygen.
        for idx in ring_set:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # only consider ring carbons
            # Look at neighbors that are not in the ring (exocyclic substituents)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring_set:
                    continue
                # Check if the neighbor is a hydroxyl group substituent
                if is_hydroxyl(nbr):
                    polar_substituent_count += 1
                # Check if the neighbor is a free -NH2 group substituent
                elif is_free_amino(nbr):
                    polar_substituent_count += 1
                    found_amino = True
        # We require that the candidate sugar ring has at least 3 polar substituents,
        # and at least one of them is an amino group.
        if found_amino and polar_substituent_count >= 3:
            return True, "Molecule contains a non-aromatic sugar-like ring (5- or 6-membered with one ring oxygen) decorated with polar substituents including a free amino group (-NH2) replacing a hydroxyl group."
    
    return False, "No candidate sugar ring with a free -NH2 substituent replacing an -OH group was detected."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "N[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O",  # 3-amino-3-deoxy-D-glucopyranose (should be True)
        "CC[C@@H]1CC(=O)c2c(C1)cc(O)c1C(=O)c3c(O[C@H]4C[C@@H](N)[C@@H](O)[C@H](C)O4)cccc3C(=O)c21"  # brasiliquinone A (should be false)
    ]
    for s in test_smiles:
        res, reason = is_amino_sugar(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")
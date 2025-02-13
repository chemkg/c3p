"""
Classifies: CHEBI:16389 ubiquinones
"""
#!/usr/bin/env python
"""
Classifies: Ubiquinones – any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone.
Ubiquinones usually carry a polyprenyl side chain at position 6.
This version inspects each 6-membered ring to verify:
    • Exactly 2 ring carbons bear a double-bonded oxygen (carbonyl group – the quinone motif)
    • Among the remaining 4 ring atoms there are exactly 2 oxygen substituents 
      (representing the 2,3-dimethoxy pattern, allowing either methoxy or hydroxy) 
      and exactly 1 methyl group (representing the 5-methyl substituent).
    • At least one prenyl (isoprene) fragment is attached to any of the ring atoms.
The program also returns the global count of prenyl (isoprenoid) units as an integer in the explanation.
"""

from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule belongs to the class of ubiquinones based on its SMILES string.
    It searches for a 6-membered benzoquinone ring with exactly two carbonyl groups and a 2,3-dimethoxy-5-methyl pattern,
    and it requires that at least one prenyl (isoprene) fragment is attached.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ubiquinone, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for an isoprene (prenyl) fragment.
    # This pattern matches one isoprene repeating unit.
    prenyl_query = Chem.MolFromSmarts("C/C=C(C)")
    # Get all prenyl matches globally (they may be part of a long polyprenyl chain)
    prenyl_matches = mol.GetSubstructMatches(prenyl_query)
    n_prenyl_global = len(prenyl_matches)
    
    # Helper functions to decide substituent type
    def is_double_bonded_oxygen(bond, nbr):
        # Returns True if the bond from a ring carbon to neighbor is a double bond and neighbor is oxygen.
        return (nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    
    def is_hydroxy(oatom):
        # Returns True if oxygen has at least one explicit hydrogen (–OH)
        # Note: sometimes H atoms are implicit; we use GetTotalNumHs() here.
        return (oatom.GetAtomicNum() == 8 and oatom.GetTotalNumHs() >= 1)
    
    def is_methoxy(oatom, from_atom_idx):
        # Check if the oxygen (that is not a double-bonded carbonyl oxygen) is linked to a methyl group.
        # We require that one of its neighbors (other than the ring atom) is a CH3 group.
        if oatom.GetAtomicNum() != 8:
            return False
        for nb in oatom.GetNeighbors():
            if nb.GetIdx() == from_atom_idx:
                continue
            # A methyl carbon should have atomic num 6, degree 1 and 3 hydrogens.
            if nb.GetAtomicNum() == 6 and nb.GetDegree() == 1 and nb.GetTotalNumHs() == 3:
                return True
        return False

    def is_methyl(catom):
        # Check if a carbon is a methyl: atomic no.6, degree=1 and 3 total hydrogens.
        return (catom.GetAtomicNum() == 6 and catom.GetDegree() == 1 and catom.GetTotalNumHs() == 3)
    
    # For each 6-membered ring, test if it matches the core pattern.
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        if len(ring) != 6:
            continue
        
        # First pass: Identify which ring atoms are carbonyl centers.
        carbonyl_indices = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # only consider carbons
            # Look at neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if is_double_bonded_oxygen(bond, nbr):
                    carbonyl_indices.add(idx)
                    break  # only mark once per ring atom
        
        # We require exactly 2 carbonyl groups on the ring.
        if len(carbonyl_indices) != 2:
            continue
        
        # For the non-carbonyl atoms in the ring, count substituents.
        oxy_substituent_count = 0   # for –OH or –OCH3 groups
        methyl_substituent_count = 0  # for direct CH3 groups
        prenyl_attached = False
        
        # To help ensure that the prenyl side chain is attached to the candidate core,
        # we record the set of non-ring neighbors of any ring atom.
        external_neighbors = set()
        
        for idx in ring:
            # Skip atoms that were already counted as having a carbonyl (their substituent is the =O).
            atom = mol.GetAtomWithIdx(idx)
            if idx in carbonyl_indices:
                continue
            # For each neighbor not in the ring:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                external_neighbors.add(nbr.GetIdx())
                # If neighbor is oxygen and not involved in a double bond (already counted as carbonyl), check substituent type.
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        continue  # skip if this were a carbonyl (should not happen here)
                    # Count if it is hydroxy or methoxy.
                    if is_hydroxy(nbr) or is_methoxy(nbr, atom.GetIdx()):
                        oxy_substituent_count += 1
                # If neighbor is carbon, check if it is a methyl group.
                elif nbr.GetAtomicNum() == 6:
                    if is_methyl(nbr):
                        methyl_substituent_count += 1
        
        # The desired pattern on the non-carbonyl ring atoms is:
        # exactly 2 oxygen substituents and exactly 1 methyl substituent.
        if oxy_substituent_count != 2 or methyl_substituent_count != 1:
            continue
        
        # Now check if at least one prenyl (isoprene unit) fragment is attached to one of the ring atoms.
        # We filter global prenyl matches to see if any has at least one atom that is directly attached
        # to a ring atom in this candidate core.
        for match in prenyl_matches:
            # For each prenyl match, check if any atom is an external neighbor of a ring atom.
            for atom_idx in match:
                if atom_idx in external_neighbors:
                    prenyl_attached = True
                    break
            if prenyl_attached:
                break
        
        # If core pattern is met and a prenyl fragment is attached, we classify as a ubiquinone.
        if prenyl_attached:
            reason = (f"Contains 2,3-dimethoxy-5-methylbenzoquinone derived core "
                      f"with {n_prenyl_global} isoprenoid unit(s) detected")
            return True, reason

    # If no candidate ring passed our screening, return false.
    return False, ("Does not contain a suitable 2,3-dimethoxy-5-methylbenzoquinone derived core "
                   "with proper substituents and a prenyl side chain")
    
# Example usage (uncomment the block below to test):
# if __name__ == "__main__":
#     tests = [
#         "S(C=1C(=O)C(OC)=C(OC)C(C1CC(=O)OC)=O)C",  # Coptirhoquinone A (expected false by our definition)
#         "COC1=C(OC)C(=O)C(CC=C(C)C)=C(C)C1=O",      # ubiquinone-1 (expected true)
#         "O=C1C(OC)=C(OC)C(=O)C(=C1CC=C(CCCC(=O)C)C)C",# Pseudoalteromone A (expected true)
#     ]
#     for s in tests:
#         result, explanation = is_ubiquinones(s)
#         print(result, explanation)
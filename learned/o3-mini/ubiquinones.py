"""
Classifies: CHEBI:16389 ubiquinones
"""
#!/usr/bin/env python
"""
Classifies: Ubiquinones – any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone.
Ubiquinones usually carry a polyprenyl side chain at position 6.
The improved algorithm checks for an isolated 6-membered ring with exactly two carbonyl groups,
two oxygen substituents (allowing -OCH3 or -OH) and one methyl group.
It also verifies that at least one isoprenoid (prenyl) fragment is attached.
"""

from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule belongs to the class of ubiquinones based on its SMILES string.
    The algorithm looks for a six‐membered ring representing the benzoquinone core that:
      – Contains exactly two carbonyl groups (C(=O)) (the quinone functionality)
      – Has exactly two oxygen substituents (which can be methoxy or hydroxy) and one methyl group 
        on the remaining ring atoms (approximating the 2,3-dimethoxy-5-methyl pattern).
      – Has at least one prenyl side chain attached to one of the ring atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a ubiquinone, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a simple SMARTS for an isoprene (prenyl) fragment.
    prenyl_smarts = "C/C=C(C)"
    prenyl_query = Chem.MolFromSmarts(prenyl_smarts)
    n_prenyl_global = len(mol.GetSubstructMatches(prenyl_query)) if prenyl_query is not None else 0
    
    # Get the ring information from the molecule.
    rings = mol.GetRingInfo().AtomRings()
    
    # Iterate over each ring looking for a candidate core.
    for ring in rings:
        # We consider only 6-membered rings.
        if len(ring) != 6:
            continue

        # Count how many carbon atoms in the ring are in a carbonyl group.
        carbonyl_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # Only carbons can be carbonyl centers
                continue
            # Look for a double-bonded oxygen neighbor
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond.GetBondTypeAsDouble() == 2.0:
                        carbonyl_indices.append(idx)
                        break
        if len(carbonyl_indices) != 2:
            continue  # need exactly 2 carbonyl groups in the ring

        # Now check the substituents on the remaining (non-carbonyl) ring atoms.
        alkoxy_count = 0  # count of O-substituents (–OCH3 or –OH)
        methyl_count = 0  # count of methyl substituents
        prenyl_attached = False
        for idx in ring:
            # Skip the carbonyl atoms.
            if idx in carbonyl_indices:
                continue
            atom = mol.GetAtomWithIdx(idx)
            # Analyze each neighbor not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check for oxygen substituent (either -OH or -OCH3)
                if nbr.GetAtomicNum() == 8:
                    # If the oxygen explicitly carries a hydrogen, accept it as hydroxyl.
                    if nbr.GetTotalNumHs() > 0:
                        alkoxy_count += 1
                    else:
                        # Alternatively, if the oxygen is bound to a carbon that appears to be a methyl group.
                        for onbr in nbr.GetNeighbors():
                            if onbr.GetIdx() == atom.GetIdx():
                                continue
                            if onbr.GetAtomicNum() == 6 and onbr.GetTotalNumHs() == 3 and onbr.GetDegree() == 1:
                                alkoxy_count += 1
                                break
                # Check for a methyl substituent (a carbon with 3 hydrogens and only one connection).
                elif nbr.GetAtomicNum() == 6:
                    if nbr.GetTotalNumHs() == 3 and nbr.GetDegree() == 1:
                        methyl_count += 1
                    # Also, check if a prenyl fragment might be attached.
                    if prenyl_query is not None:
                        # If the neighbor (or the connection) is part of a prenyl fragment, mark it.
                        matches = mol.GetSubstructMatches(prenyl_query)
                        for match in matches:
                            if atom.GetIdx() in match or nbr.GetIdx() in match:
                                prenyl_attached = True
                                break
            
        # We require at least two oxygen substituents and one methyl on the non-carbonyl positions,
        # plus that one of the ring atoms has a prenyl side chain attached.
        if alkoxy_count >= 2 and methyl_count >= 1 and prenyl_attached:
            reason = (f"Contains 2,3-dimethoxy-5-methylbenzoquinone derived core "
                      f"with {n_prenyl_global} isoprenoid unit(s) detected")
            return True, reason

    return False, ("Does not contain a suitable 2,3-dimethoxy-5-methylbenzoquinone derived core "
                   "with proper substituents and a prenyl side chain")
    
# Example usage (uncomment below to test):
# if __name__ == "__main__":
#     test_smiles = "COC1=C(OC)C(=O)C(CC=C(C)C)=C(C)C1=O"  # ubiquinone-1
#     result, explanation = is_ubiquinones(test_smiles)
#     print(result, explanation)
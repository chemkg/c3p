"""
Classifies: CHEBI:18000 aralkylamine
"""
#!/usr/bin/env python3
"""
Classifies: aralkylamine
Definition: An alkylamine in which the alkyl substituent is itself substituted by an aromatic group.
That is, a non‐aromatic (aliphatic) nitrogen atom that is connected (via at least one, but not directly)
to an alkyl (carbon only) chain that eventually bears an aromatic ring.
"""

from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine having an alkyl substituent that bears an aromatic group.
    More specifically, the amine nitrogen should not be directly bonded to an aromatic ring 
    (that would be an arylamine) and should not be part of an amide.
    Instead, it must have at least one bond (via a chain of non‐aromatic carbons) leading to an aromatic atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple containing True/False and a reason.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: DFS that starts from a carbon (found as an alkyl substituent of N)
    # and only traverses non‐aromatic carbons.
    # It returns True if an aromatic atom is reached within max_depth bonds.
    def dfs_aralkyl(atom, distance, max_distance, visited):
        # We only want an aromatic hit if it is not directly attached to the N.
        # (i.e. require at least one intervening carbon bond).
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            # If the neighbor is aromatic and we are at least one bond away from the starting carbon,
            # then we have found an aromatic substituent.
            if nbr.GetIsAromatic() and (distance + 1 >= 1):
                return True
            # Only continue the search on non‐aromatic carbons.
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()):
                # Only traverse further if we still have depth remaining.
                if distance + 1 < max_distance:
                    visited.add(nbr.GetIdx())
                    if dfs_aralkyl(nbr, distance + 1, max_distance, visited):
                        return True
        return False

    # Iterate over all atoms looking for candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        # Look for nitrogen atoms that are non‐aromatic (i.e. typical alkylamine N)
        if atom.GetAtomicNum() == 7 and (not atom.GetIsAromatic()):
            # Exclude N atoms that are part of an amide function.
            # That is, if the N is bonded to a carbon which is double‐bonded to an oxygen,
            # then skip this nitrogen.
            is_amide = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:  # candidate carbon
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    for nbr2 in nbr.GetNeighbors():
                        # Check if this carbon is bonded to an oxygen with a double bond.
                        if nbr2.GetAtomicNum() == 8:
                            bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                            if bond2 is not None and bond2.GetBondType() == Chem.BondType.DOUBLE:
                                is_amide = True
                                break
                    if is_amide:
                        break
            if is_amide:
                # Skip this nitrogen since it is in an amide (or similar) context.
                continue

            # For each neighbor of the candidate N, check if it is a carbon chain that leads to an aromatic group.
            for nbr in atom.GetNeighbors():
                # We insist that the chain begins with a carbon (not an aromatic one)
                # so that we do not confuse a direct aromatic bond (arylamine) as a positive.
                if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()):
                    # We use a DFS from this carbon. Set a max allowed chain length (here 3 bonds).
                    if dfs_aralkyl(nbr, distance=0, max_distance=3, visited={nbr.GetIdx()}):
                        return True, "Contains an aliphatic amine nitrogen with an alkyl substituent leading to an aromatic group"
    
    return False, "No aralkyl substituent detected on an alkylamine nitrogen"

# For testing purposes, you may run the following examples:
if __name__ == "__main__":
    # Benzylamine is the simplest example: N attached to CH2 then aromatic ring ("NCc1ccccc1")
    test_smiles = "NCc1ccccc1"
    result, reason = is_aralkylamine(test_smiles)
    print("SMILES:", test_smiles)
    print("Aralkylamine?", result)
    print("Reason:", reason)
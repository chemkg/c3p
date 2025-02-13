"""
Classifies: CHEBI:18000 aralkylamine
"""
#!/usr/bin/env python3
"""
Classifies: aralkylamine
Definition: An alkylamine in which at least one alkyl substituent (i.e. a chain of only sp3 carbons)
leads (indirectly) to an aromatic ring. In other words, an aliphatic (non‐aromatic) amine N that is 
connected via at least one “alkyl” (saturated carbon) bond to an aromatic ring.
This excludes cases where the nitrogen is directly bonded to an aromatic ring (arylamines)
or when it is in an amide context.
"""

from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine having an alkyl substituent (i.e. only sp3 carbons)
    that eventually carries an aromatic ring. The amine nitrogen should not be directly bonded to an 
    aromatic ring, and should not be part of an amide linkage.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the boolean is True if the molecule contains at least one aralkylamine,
        and the string is a reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude non-relevant cases: e.g. if the molecule has no nitrogen:
    if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "No nitrogen atoms present"

    # Recursive search: from a starting sp3 non‐aromatic carbon,
    # traverse only along sp3 carbon atoms (i.e. pure alkyl chain) to look for an aromatic atom.
    # We allow a chain length up to max_depth bonds.
    def dfs_chain(atom, distance, max_depth, visited):
        # For each neighbor, we only traverse if it is a pure (non‐aromatic) carbon.
        for nbr in atom.GetNeighbors():
            # If this neighbor is aromatic, then we have reached an aromatic substituent.
            # We require that at least one bond has been traversed from the N (distance>=0 means starting carbon already)
            if nbr.GetIsAromatic() and (distance + 1 >= 1):
                return True
            # Only continue if the neighbor is carbon, not aromatic, and sp3 hybridized.
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()):
                # Check hybridization: we want only typical sp3 alkyl carbons.
                if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                # Avoid going backwards.
                if nbr.GetIdx() in visited:
                    continue
                if distance + 1 < max_depth:
                    visited.add(nbr.GetIdx())
                    if dfs_chain(nbr, distance + 1, max_depth, visited):
                        return True
        return False

    # Now iterate over all candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and (not atom.GetIsAromatic()):
            # Exclude N atoms that are part of an amide (or similar)
            is_amide = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:  # candidate carbon attached to N
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is None:
                        continue
                    # Check if that carbon is double-bonded to an oxygen (typical of amides)
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 8:
                            bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                            if bond2 is not None and bond2.GetBondType() == Chem.BondType.DOUBLE:
                                is_amide = True
                                break
                    if is_amide:
                        break
            if is_amide:
                continue  # Skip candidate nitrogen in an amide context

            # For each neighbor of candidate N: we require that it is a carbon and is not aromatic.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()):
                    # Also, check that the neighbor is sp3. This ensures that the chain is an alkyl chain.
                    if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                        continue
                    # We want to ensure that N is not directly bonded to an aromatic ring.
                    found = dfs_chain(nbr, distance=0, max_depth=6, visited={nbr.GetIdx()})
                    if found:
                        return True, "Contains an aliphatic amine nitrogen with an alkyl substituent leading to an aromatic ring"
    return False, "No appropriate aralkyl substituent found on an aliphatic amine nitrogen"

# For testing purposes, one can run some examples.
if __name__ == "__main__":
    # Simple benzylamine example: N attached to CH2 then aromatic ring.
    test_smiles = "NCc1ccccc1"
    result, reason = is_aralkylamine(test_smiles)
    print("SMILES:", test_smiles)
    print("Aralkylamine?", result)
    print("Reason:", reason)
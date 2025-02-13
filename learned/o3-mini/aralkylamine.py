"""
Classifies: CHEBI:18000 aralkylamine
"""
#!/usr/bin/env python3
"""
Classifies: aralkylamine
Definition: An alkylamine in which at least one alkyl substituent (i.e. a chain of only sp3 carbons)
leads (indirectly) to an aromatic ring. That is, the molecule contains a non‐aromatic (aliphatic) amine N,
which is not directly bound to an aromatic ring and not part of an amide, and which is connected via at least
one saturated carbon (alkyl) bond to an aromatic carbon.
"""

from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    An aralkylamine is defined as an alkylamine having at least one alkyl (saturated, sp3 hybridized)
    substituent that eventually leads to an aromatic ring. In this context the substitution must be via
    at least one saturated carbon. Additionally, the nitrogen itself must not be directly bonded to an
    aromatic ring (which would be an arylamine) and should not be involved in an amide linkage.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the boolean is True if the molecule contains at least one aralkylamine,
        and the string provides a reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If no nitrogen is present, this is not an aralkylamine.
    if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "No nitrogen atoms present"
    
    # Define a helper to check if a nitrogen is attached to a carbon that is double-bonded to oxygen
    # (i.e. likely part of an amide or similar).
    def in_amide_context(nitrogen):
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # candidate carbon
                bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # Check if the carbon has a double bond with oxygen (amide C=O)
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 8:
                        bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bond2 is not None and bond2.GetBondType() == Chem.BondType.DOUBLE:
                            return True
        return False

    # DFS function:
    # From a given carbon atom in the alkyl chain (already guaranteed to be sp3 and non-aromatic),
    # we traverse along neighbors that are also sp3 carbons and non-aromatic.
    # If we eventually encounter an aromatic carbon (which means the chain leads to an aromatic ring),
    # we return True.
    def dfs_chain(atom, current_depth, max_depth, visited):
        # Traverse neighbors of the current atom
        for nbr in atom.GetNeighbors():
            # If neighbor is a carbon atom
            if nbr.GetAtomicNum() == 6:
                # If the neighbor is aromatic, we have reached an aromatic substituent.
                if nbr.GetIsAromatic():
                    # We allow a match if we have traversed at least one saturated alkyl carbon.
                    if current_depth >= 1:
                        return True
                # Traverse if neighbor is saturated sp3 and non-aromatic.
                elif nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and not nbr.GetIsAromatic():
                    idx = nbr.GetIdx()
                    if idx in visited:
                        continue
                    # Extend the chain if not too long.
                    if current_depth < max_depth:
                        visited.add(idx)
                        if dfs_chain(nbr, current_depth + 1, max_depth, visited):
                            return True
        return False
    
    # Iterate over all candidate nitrogen atoms that are non-aromatic.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic():
            # Skip if this nitrogen appears to be in an amide context.
            if in_amide_context(atom):
                continue
            
            # Check all neighbors of this nitrogen.
            for nbr in atom.GetNeighbors():
                # We require that the substituent starts with a saturated, non-aromatic carbon.
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() \
                   and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                    # Start DFS search from this carbon with initial depth = 1 (one bond from N).
                    visited = set([nbr.GetIdx()])
                    if dfs_chain(nbr, current_depth=1, max_depth=6, visited=visited):
                        return True, "Contains an aliphatic amine nitrogen with an alkyl substituent leading to an aromatic ring"
    return False, "No appropriate aralkyl substituent found on an aliphatic amine nitrogen"


# For testing purposes.
if __name__ == "__main__":
    # Example: benzylamine is an aralkylamine because N-CH2-c1ccccc1 (here the chain goes CH2 then aromatic)
    test_smiles = "NCc1ccccc1"
    result, reason = is_aralkylamine(test_smiles)
    print("SMILES:", test_smiles)
    print("Aralkylamine?", result)
    print("Reason:", reason)
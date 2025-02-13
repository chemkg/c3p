"""
Classifies: CHEBI:24026 fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: Fatty alcohol
Definition: An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.
We check that:
 1. The molecule has at least one hydroxyl group (–OH) attached to a non‐aromatic (aliphatic) carbon.
 2. Starting from that carbon, we determine the length of the contiguous chain of non‐aromatic carbon atoms.
    If the longest chain has at least 3 carbon atoms, we classify it as a fatty alcohol.
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is defined as an aliphatic alcohol having a chain of at least 3 carbon atoms.
    (The chain can be saturated, unsaturated, and may be branched.)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Identify hydroxyl groups: Oxygen with one hydrogen neighbor.
    # We then check that the -OH is attached to a carbon and that that carbon is not aromatic.
    hydroxyl_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen
            # Check that it is in an -OH group (has 1 hydrogen neighbor)
            hydrogens = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
            if len(hydrogens) != 1:
                continue
            # Now check its neighbor: Must be a carbon.
            carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if not carbons:
                continue
            # For our purposes we require the carbon is aliphatic (not aromatic)
            for c in carbons:
                if not c.GetIsAromatic():
                    hydroxyl_atoms.append(c)
    if not hydroxyl_atoms:
        return False, "No aliphatic hydroxyl group found"

    # Define a helper function to compute the longest simple chain of aliphatic (non-aromatic) carbons.
    def dfs(atom, visited):
        """
        Depth-first search returning the length of the longest chain starting at the given carbon atom.
        Only considers neighbors that are carbon and non-aromatic.
        """
        max_length = 1  # count this atom
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (nbr.GetIdx() not in visited):
                length = 1 + dfs(nbr, visited.copy())  # use a copy so that each branch is independent
                if length > max_length:
                    max_length = length
        return max_length

    # Check each candidate aliphatic -OH attached carbon: get the longest aliphatic (non-aromatic C) chain.
    for c_atom in hydroxyl_atoms:
        chain_length = dfs(c_atom, set())
        # According to the definition, the chain length must be at least 3 carbons.
        if chain_length >= 3:
            reason = f"Found aliphatic chain of {chain_length} carbon atoms attached to the hydroxyl group."
            return True, reason

    return False, "No contiguous aliphatic carbon chain of minimum length (3) found attached to hydroxyl group."

# For testing purposes:
if __name__ == "__main__":
    # List of example SMILES strings that are supposed to be fatty alcohols.
    examples = [
        "CCCCCCCCCCCCCCC(O)CCCC",  # nonadecan-5-ol
        "O(C(CCCCCCCCCC(O)C)CC1=CC(O)=CC(O)=C1)C(=O)C",  # 1-(3,5-Dihydroxyphenyl)-12-hydroxytridecan-2-yl acetate
        "CCCCCCCCCCCCCC(O)CCCC",  # octadecan-5-ol
        "OC(CCCCC)CC(=O)CCCCCCCCCCCCC",  # 6-Hydroxy-8-heneicosanone
        "CCCCCCCCCCCCCCC(O)CCCCCCCCCC",  # pentacosan-11-ol
        "OC(C(O)CO)C(O)C#CC#CC#CC#CC",  # (2S,3S,4S)-5,7,9,11-Tridecatetrayne-1,2,3,4-tetrol
        "CCCCCCC(O)CCC",  # undecan-4-ol (simplified)
        "CCCCCCCCCCC(O)CCCCC",  # hexadecan-6-ol
        "CCCCCCC(O)CCCCCCCO",  # 1,8-tetradecanediol
        "CCCCCCCCCCCC(O)CCCCCC",  # nonadecan-7-ol
        "CC(C)CCCCCCCCCCCCO",  # 13-methyltetradecan-1-ol
        "CCCCCCCCCC(O)CCCCCCC",  # hexadecan-8-ol
        "OCC/C=C\\CCCCCCCC/C=C\\CCCC",  # 3Z,13Z-octadecadien-1-ol
        "CCCCCCCCCCCCCCCCCC(O)CC",  # henicosan-3-ol
        "C(CCCCCC)O",  # heptan-1-ol
        "C(/C=C/CC/C=C\\CC)O",  # (2E,6Z)-nona-2,6-dien-1-ol
        "CCCCCCCCCCCCCCCC(O)CCCCCCC",  # tricosan-8-ol
        "O[C@H](CCCCCCCCCCCCCCCC)COC[C@@H](O)CO",  # 1-O-(2R-hydroxy-octadecyl)-sn-glycerol
        "O=C(C1=C(O)C=C(O)C=C1O)CCC/C=C\\C/C=C\\C\\C=C/C=C/[C@@H](O)C/C=C\\CC",  # 15-Hydroxy-1-(2,4,6-trihydroxyphenyl)-5Z,8Z,11Z,13E,17Z-eicosapentaen-1-one
        "CCCCCCCCCCCC(O)CCCCCCCCCC",  # docosan-11-ol
        "OCCCCCCCCCCCCC#CCC=C",  # 17-OCTADECEN-14-YN-1-OL
        "CCC(C)CCCCCCCCCCCCCO",  # 14-methylhexadecan-1-ol
        "CCCCCCCCCCCCC(O)CCCCCCCC",  # octadecan-9-ol
        "CCCCCCCCCC(O)CCCCCC",  # hexadecan-7-ol
        "CCCCCCCCC(O)",  # nonan-1-ol
        "CCCC(O)CCCCCCCCCCCC",  # 1,7-decanediol (simplified)
    ]
    for s in examples:
        result, reason = is_fatty_alcohol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")
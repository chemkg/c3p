"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine

Definition:
An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl groups are located 
at positions 1 and 2 respectively. Such molecules comprise three parts:
    • An O-alkyl (ether) chain (position 1)
    • An O-acyl (ester) chain (position 2)
    • A phosphocholine headgroup (attached via phosphate at position 3)

This program uses RDKit substructure searches and additional heuristics to distinguish this lipid type 
from diacyl phosphocholines.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_alkyl_chain(atom, visited=None):
    """
    Recursively find the longest path (in number of carbon atoms) in an acyclic, aliphatic chain
    starting from the given atom. Only traverses carbons that are aliphatic (atomic num 6, sp3) 
    and not involved in carbonyl (as determined by whether a double bond to oxygen exists).
    """
    if visited is None:
        visited = set()
    visited.add(atom.GetIdx())
    max_length = 0
    # Look at neighbors that are carbons
    for nbr in atom.GetNeighbors():
        # Only traverse carbon atoms that are sp3 (and not in rings) and are not yet visited.
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
            # Exclude carbons that are carbonyl: if any bond from nbr is a double bond to O
            is_carbonyl = False
            for bond in nbr.GetBonds():
                if bond.GetBondTypeAsDouble() >= 2:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8:
                        is_carbonyl = True
                        break
            if is_carbonyl:
                continue
            branch_length = 1 + longest_alkyl_chain(nbr, visited.copy())
            if branch_length > max_length:
                max_length = branch_length
    return max_length

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Improved heuristic logic:
      1. Confirm the presence of the phosphocholine headgroup. We look for an O-phosphate 
         substructure that is directly connected to "CC[N+](C)(C)C".
      2. Count the number of ester groups represented by the SMARTS "OC(=O)". The target lipids 
         should have exactly one acyl (ester) chain.
      3. Identify an O-alkyl (ether) group that is not bound to phosphorus and not part of an ester.
         For candidate oxygens, check that at least one of the neighboring carbons leads to a long (>=8)
         uninterrupted alkyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule passes the tests and is deemed a member of this lipid class.
        str: A reason string for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for phosphocholine headgroup.
    # Pattern: an oxygen attached to phosphate which is further connected to "CC[N+](C)(C)C"
    phos_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phos_pattern):
        return False, "Missing phosphocholine headgroup (pattern OP(=O)([O-])OCC[N+](C)(C)C not found)"

    # 2. Count ester (acyl) groups using the SMARTS pattern "OC(=O)"
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "Missing ester (acyl) group (pattern OC(=O) not found)"
    if len(ester_matches) != 1:
        return False, f"Expected exactly one acyl ester group but found {len(ester_matches)}; might be diacyl."

    # 3. Look for an O-alkyl (ether) linkage.
    # We iterate over all oxygen atoms:
    ether_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue

        # Skip if oxygen is part of the phosphate headgroup: check if any neighbor is phosphorus.
        if any(neigh.GetAtomicNum() == 15 for neigh in atom.GetNeighbors()):
            continue

        # Skip if oxygen is part of an ester:
        # Check if one neighbor is a carbon that has a double bond to oxygen (carbonyl).
        part_of_ester = False
        neighbor_carbons = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                neighbor_carbons.append(nbr)
                for bond in mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetOtherAtoms([atom])[1:]:
                    # Not the safest method; instead, check bonds at neighbor:
                    for b in nbr.GetBonds():
                        if b.GetBondTypeAsDouble() >= 2:
                            other = b.GetOtherAtom(nbr)
                            if other.GetAtomicNum() == 8:
                                part_of_ester = True
                                break
                    if part_of_ester:
                        break
            if part_of_ester:
                break
        if part_of_ester:
            continue  # skip oxygens in an ester pattern

        # Now, for each remaining oxygen, check that at least one carbon neighbor leads to a long alkyl chain
        for nbr in neighbor_carbons:
            # Calculate the length of the longest linear alkyl chain starting at this carbon.
            chain_length = longest_alkyl_chain(nbr)
            if chain_length >= 8:
                ether_found = True
                break
        if ether_found:
            break

    if not ether_found:
        return False, "Missing or too short O-alkyl (ether) chain; no oxygen found leading to alkyl chain of sufficient length."

    return True, ("Contains phosphocholine headgroup, exactly one acyl ester group, and one O-alkyl (ether) chain "
                   "with a long carbon chain. This is consistent with 2-acyl-1-alkyl-sn-glycero-3-phosphocholine.")

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCC"
    result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(test_smiles)
    print(f"Test result: {result}\nReason: {reason}")
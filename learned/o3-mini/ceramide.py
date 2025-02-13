"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)
Definition:
  Ceramides are sphingoid bases with an amide-linked fatty acid.
  The fatty acid is typically saturated or monounsaturated with chain lengths
  from 14 to 26 carbon atoms; the sphingoid base usually carries at least one hydroxyl group.
  
This implementation:
  - Iterates through all C(=O)-N bonds (amide bonds) in the molecule.
  - For each amide bond, it identifies the acyl chain on the carbonyl carbon,
    measures its maximum contiguous carbon chain length (heuristically),
    and requires this to be between 14 and 26.
  - For the sphingoid side (the substituent on the nitrogen away from the carbonyl),
    it verifies that at least one oxygen is present (to mimic a common hydroxyl group on C2).
    
Note:
  This is a heuristic method and may not capture all ceramide structures.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    A ceramide is defined as an N-acyl-sphingoid base:
      - It must contain an amide bond (C(=O)N) linking a fatty acid (14-26 carbons) 
        to a sphingoid base.
      - The sphingoid base side should have at least one oxygen (to indicate a hydroxyl).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a ceramide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string to get a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: Recursively compute the maximum length of a contiguous chain of carbons.
    # This function traverses only through carbon atoms.
    def get_max_chain_length(atom, coming_from_idx, visited):
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == coming_from_idx:
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                chain_length = 1 + get_max_chain_length(nbr, atom.GetIdx(), visited.copy())
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # Iterate through all bonds to find potential amide bonds.
    for bond in mol.GetBonds():
        # Check if this bond connects carbon and nitrogen.
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7) or (a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6):
            # Identify the carbon and nitrogen atoms.
            if a1.GetAtomicNum() == 6:
                carbonyl_atom = a1
                nitrogen_atom = a2
            else:
                carbonyl_atom = a2
                nitrogen_atom = a1
            
            # Verify that the carbonyl_atom bears a double-bonded oxygen (i.e. is a carbonyl)
            has_double_oxygen = False
            for nbr in carbonyl_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond_to_oxygen = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), nbr.GetIdx())
                    if bond_to_oxygen.GetBondType() == Chem.BondType.DOUBLE:
                        has_double_oxygen = True
                        break
            if not has_double_oxygen:
                continue  # Not a carbonyl carbon
            
            # Identify the acyl chain: a neighbor of the carbonyl carbon that is not the nitrogen
            acyl_candidate = None
            for nbr in carbonyl_atom.GetNeighbors():
                if nbr.GetIdx() == nitrogen_atom.GetIdx():
                    continue
                # Skip the oxygen from the carbonyl
                if nbr.GetAtomicNum() == 8:
                    continue
                if nbr.GetAtomicNum() == 6:
                    acyl_candidate = nbr
                    break
            if acyl_candidate is None:
                continue  # Could not find an acyl chain attachment
            
            # Calculate the acyl chain length using the helper function.
            acyl_chain_length = get_max_chain_length(acyl_candidate, carbonyl_atom.GetIdx(), set())
            if not (14 <= acyl_chain_length <= 26):
                # Fatty acid chain length out of ceramide range.
                continue
            
            # Identify the sphingoid side: on nitrogen, choose a neighbor that is not the carbonyl carbon.
            sphingo_candidate = None
            for nbr in nitrogen_atom.GetNeighbors():
                if nbr.GetIdx() == carbonyl_atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    sphingo_candidate = nbr
                    break
            if sphingo_candidate is None:
                continue  # No candidate sphingoid base found

            # Check if the sphingoid candidate contains at least one oxygen (as a proxy for a hydroxyl group)
            sphingo_has_OH = False
            for nbr in sphingo_candidate.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    sphingo_has_OH = True
                    break
            if not sphingo_has_OH:
                continue  # Sphingoid part lacks oxygen/hydroxyl group
            
            # If we reach here, we have a matching amide bond with the expected acyl chain and sphingoid OH.
            return True, f"Ceramide detected; acyl chain length = {acyl_chain_length} carbons."
    
    return False, "No amide bond with both an acyl chain (14-26 C) and a sphingoid side possessing a hydroxyl was detected."

# (Optional) Basic testing code.
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)C(O)CCCC",  # Example 1
        "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCC(C)C",  # Example 2
        "CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCC"  # Example 3
    ]
    for smi in test_smiles:
        result, reason = is_ceramide(smi)
        print("SMILES:", smi)
        print("is_ceramide:", result)
        print("Reason:", reason)
        print()
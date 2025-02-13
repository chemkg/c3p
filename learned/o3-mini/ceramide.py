"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)
Definition:
  Ceramides are sphingoid bases with an amide-linked fatty acid. The fatty acid is typically saturated or 
  monounsaturated with 14–26 carbon atoms and the sphingoid base usually has a hydroxyl group (often at C2).

This program uses an approximate substructure pattern to identify a ceramide:
  - It looks for the substructure "C(=O)N[C](CO)" where the amide carbonyl is linked to 
    a nitrogen that is in turn connected to a carbon with an attached –OH (the "CO" fragment).
  - It then inspects the acyl chain (the carbon connected to the carbonyl carbon, not part of the amide)
    and measures its length. The chain length must be between 14 and 26 carbons.
  - Also, it verifies that the sphingoid side (the carbon attached to the amide nitrogen on the head-group)
    carries at least one oxygen atom (as a proxy for a hydroxyl group).
    
Note: This is a heuristic method and may not cover all edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    A ceramide is defined as an N-acyl-sphingoid base:
      - It must contain an amide bond (C(=O)N) linking a fatty acid to a sphingoid base.
      - The fatty acid chain (attached to the carbonyl carbon) should be between 14 and 26 carbons.
      - The sphingoid base side should exhibit at least one hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a ceramide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a substructure pattern having an amide and a nearby hydroxymethyl (CO) group.
    # This pattern is expected to capture many ceramide cores.
    # Pattern: C(=O)N[C](CO)
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C](CO)")
    matches = mol.GetSubstructMatches(ceramide_pattern)
    if not matches:
        return False, "No ceramide core substructure (C(=O)N[C](CO)) detected"
    
    # Helper function: recursively determine the maximum length of a contiguous chain of carbons.
    def get_max_chain_length(atom, coming_from_idx, visited):
        # Only traverse if the atom is a carbon.
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1  # count this atom
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == coming_from_idx:
                continue
            # Continue only if the neighbor is carbon and not visited (to avoid cycles)
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # Recursive call passing current atom's index as coming from.
                chain_length = 1 + get_max_chain_length(nbr, atom.GetIdx(), visited.copy())
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # Loop over each substructure match.
    for match in matches:
        # According to our SMARTS pattern, the match tuple is:
        # match[0] = carbonyl carbon "C(=O)", match[1] = amide nitrogen "N", match[2] = sphingoid carbon with CO.
        carbonyl_idx, amideN_idx, sphingo_idx = match[0], match[1], match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        amideN_atom = mol.GetAtomWithIdx(amideN_idx)
        sphingo_atom = mol.GetAtomWithIdx(sphingo_idx)
        
        # Identify the acyl chain starting from the carbonyl carbon.
        # The carbonyl carbon in an amide is bonded to:
        #   - an oxygen (from the C=O)
        #   - the nitrogen (amide bond)
        #   - and one additional substituent that should be the start of the fatty acid chain.
        acyl_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip the oxygen (usually symbol "O") and the amide nitrogen.
            if nbr.GetIdx() == amideN_idx:
                continue
            if nbr.GetAtomicNum() == 8:
                continue
            # Expecting a carbon here; pick the first one.
            if nbr.GetAtomicNum() == 6:  
                acyl_start = nbr
                break
        
        if acyl_start is None:
            continue  # skip this match if no acyl chain was found.
        
        # Count the length of the acyl chain.
        acyl_chain_length = get_max_chain_length(acyl_start, carbonyl_idx, set())
        if not (14 <= acyl_chain_length <= 26):
            continue  # chain length not in the expected range.
        
        # Now verify the sphingoid side: Check that the nitrogen has another substituent (the sphingo head)
        # which should have at least one oxygen attached (as a proxy for a hydroxyl group).
        sphingo_candidate = None
        for nbr in amideN_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue
            # We expect this neighbor to be a carbon: the start of the sphingoid base.
            if nbr.GetAtomicNum() == 6:
                sphingo_candidate = nbr
                break
        if sphingo_candidate is None:
            continue

        # Check if the sphingoid candidate has at least one oxygen neighbor.
        has_OH = False
        for nbr in sphingo_candidate.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                has_OH = True
                break
        if not has_OH:
            continue  # No hydroxyl on the sphingoid base.
        
        # If we reach here, we have found a ceramide core with appropriate acyl chain length and a hydroxyl.
        return True, f"Ceramide detected with acyl chain length {acyl_chain_length} carbons"
    
    # If none of the matches fulfill our criteria, return False.
    return False, "No substructure had both an amide with appropriate acyl chain length and a sphingoid hydroxyl group"

# (Optional) basic testing
if __name__ == "__main__":
    # A list of sample ceramide SMILES (some from the provided examples)
    test_smiles = [
        "CCCCCCCCCCCCCCCC[C@@H](O)[C@H](O)\\C=C\\CCCCCCCCCC(C)C",  # N-heptadecanoyl-15-methylhexadecasphing-4-enine
        "CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCC"  # N-octacosanoyl-15-methylhexadecasphinganine
    ]
    for smi in test_smiles:
        res, reason = is_ceramide(smi)
        print(f"SMILES: {smi}\n  is_ceramide: {res}\n  Reason: {reason}\n")
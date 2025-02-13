"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax Ester
Definition: A fatty acid ester resulting from the condensation of the carboxy group 
of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.
This implementation ensures that exactly one ester group is present, breaks the molecule 
at the ester linkage, and thereafter verifies that the two resulting fragments (fatty acid and 
fatty alcohol) are acyclic, predominantly carbon in composition, and long enough.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a simple wax ester based on its SMILES string.
    A wax ester is defined as a fatty acid ester (a fatty acid and a fatty alcohol
    connected via a single ester linkage).
    
    The molecule is first checked for exactly one ester group. Then the bond between the 
    ester oxygen and its adjacent carbon (on the fatty alcohol side) is broken so that the 
    two fragments can be isolated and examined for acyclic, aliphatic character.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a wax ester, False otherwise.
        str: Explanation message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an ester pattern for a simple ester: carbonyl carbon connected to an oxygen with no hydrogen.
    ester_smarts = "[#6X3](=O)[O;H0]"  
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found that meets the wax ester pattern"
    if len(ester_matches) > 1:
        return False, "More than one ester group found; not a simple wax ester"
    
    # Unpack the match: by convention our SMARTS gives indices:
    # match[0] = carbonyl carbon, match[1]= carbonyl oxygen, match[2]=ester oxygen.
    match = ester_matches[0]
    acid_c_idx = match[0]
    ester_oxygen_idx = match[2]

    # From the ester oxygen, we need the fatty alcohol side (the non-carbonyl neighbor)
    oxy_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
    alcohol_neighbor = None
    for nbr in oxy_atom.GetNeighbors():
        if nbr.GetIdx() != acid_c_idx:
            alcohol_neighbor = nbr
            break
    if alcohol_neighbor is None:
        return False, "Could not determine fatty alcohol side of the ester"

    alcohol_c_idx = alcohol_neighbor.GetIdx()

    # Ensure that neither side of the ester is in a ring.
    if mol.GetAtomWithIdx(acid_c_idx).IsInRing() or oxy_atom.IsInRing():
        return False, "Ester group is in a ring; not a simple wax ester"
        
    # Find the bond between the ester oxygen and the alcohol side.
    bond = mol.GetBondBetweenAtoms(ester_oxygen_idx, alcohol_c_idx)
    if bond is None:
        return False, "Unable to find the bond corresponding to the ester linkage"

    bond_idx = bond.GetIdx()

    # Fragment the molecule at the ester bond.
    # This will leave dummy atoms ([*]) at the breakpoints.
    fragmented = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    frags = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
    if len(frags) != 2:
        # If fragmentation did not yield two fragments, we cannot easily assign parts.
        return False, "Fragmentation did not yield two clear parts; ambiguous structure"

    # Identify which fragment is the fatty acid side and which is the fatty alcohol.
    # We do this by checking which fragment contains acid_c_idx and which contains alcohol_c_idx.
    # Note: after fragmentation, atom indices are changed. So we match by substructure:
    acid_frag = None
    alcohol_frag = None
    for frag in frags:
        atom_ids = [atom.GetAtomMapNum() for atom in frag.GetAtoms() if atom.HasProp("molAtomMapNumber")]
        # Instead, we check if the fragment contains a substructure similar to the ester carbonyl (C=O)
        # or if it contains a dummy atom at one end.
        # A simple heuristic is to check if the fragment has an oxygen double-bonded to a carbon.
        has_carbonyl = any(a.GetSymbol()=='O' for a in frag.GetAtoms() for b in a.GetBonds() 
                           if b.GetEndAtom(a).GetSymbol()=='C' and b.GetBondTypeAsDouble() == 2.0)
        # Next, check for dummy atoms:
        has_dummy = any(a.GetAtomicNum()==0 for a in frag.GetAtoms())
        if has_carbonyl:
            acid_frag = frag
        elif has_dummy:
            alcohol_frag = frag
    # If assignment by heuristic fails, simply assign by fragment size.
    if acid_frag is None or alcohol_frag is None:
        if frags[0].GetNumAtoms() >= frags[1].GetNumAtoms():
            acid_frag, alcohol_frag = frags[0], frags[1]
        else:
            acid_frag, alcohol_frag = frags[1], frags[0]

    # Define helper function to check if a fragment is a simple, acyclic aliphatic chain.
    def validate_chain(fragment, allowed=set([6]), extra_allowed=set(), min_carbons=6):
        # allowed: atomic numbers that are normally allowed (e.g. carbon, here 6)
        # extra_allowed: additional elements permitted (for acid side an oxygen might be allowed)
        # Exclude dummy atoms (atomic num 0).
        # Also reject if the fragment contains any rings.
        if rdmolops.GetSSSR(fragment) > 0:
            return False, "Fragment is cyclic"
        carbon_count = 0
        for atom in fragment.GetAtoms():
            anum = atom.GetAtomicNum()
            if anum == 0:
                continue  # ignore dummy atoms
            if anum == 6:
                carbon_count += 1
            else:
                if anum not in extra_allowed:
                    return False, f"Fragment contains disallowed heteroatom: {atom.GetSymbol()}"
        if carbon_count < min_carbons:
            return False, f"Fragment too short; found only {carbon_count} carbon(s)"
        return True, carbon_count

    # For the fatty acid (acyl) fragment, we allow carbon (6) and one oxygen (from the carbonyl).
    valid_acid, acid_count_or_msg = validate_chain(acid_frag, allowed=set([6]), extra_allowed=set([8]), min_carbons=6)
    if not valid_acid:
        return False, "Fatty acid fragment rejected: " + acid_count_or_msg
    # For the fatty alcohol fragment, we expect only carbons.
    valid_alcohol, alcohol_count_or_msg = validate_chain(alcohol_frag, allowed=set([6]), extra_allowed=set(), min_carbons=6)
    if not valid_alcohol:
        return False, "Fatty alcohol fragment rejected: " + alcohol_count_or_msg

    message = ("CORRECT Wax ester detected: fatty acid chain with {} C(s) and fatty alcohol chain with {} C(s) "
               "connected via a single ester linkage.").format(acid_count_or_msg, alcohol_count_or_msg)
    return True, message

# Example usage:
if __name__ == "__main__":
    # Test one of the provided wax ester examples:
    test_smiles = "O(CCCCCCCCCCCCCCCCCCCC(C)C)C(CCCCCCC/C=C\\CCCCCCCC)=O"  # 1-O-20-methylhenicosyl oleate
    result, reason = is_wax_ester(test_smiles)
    print(result, reason)
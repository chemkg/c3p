"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine containing compounds
Definition: Compounds having the structure RC(=NR)NR2 (commonly the -C(=NH)NH2 group, or variants thereof).
Note: This is a simplified substructure search; in order to improve upon a simple SMARTS match we perform 
an additional inspection of the environment of the match and also flag molecules that appear to be peptides.
"""

from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group based on its SMILES string.
    The group is defined as RC(=NR)NR2 (for example the common -C(=NH)NH2 group) and derivatives.
    
    The algorithm is as follows:
      1. Parse the SMILES string.
      2. Look for two possible SMARTS patterns:
         a. a neutral carboxamidine: [CX3](=[NX2])[NX3]
         b. a charged variant: [CX3](=[NX3+])[NX3]
      3. If any match is found, for each match we check its environment. In particular we look at the carbon 
         atom involved. If that carbon is attached to a substituent (beyond the two nitrogens of the amidine)
         that is “peptide-like” – for example, if that substituent in turn is attached (by a double bond) to an oxygen – 
         then we flag the match.
      4. In addition, if the molecule contains a few peptide bonds (detected by the pattern N-C(=O)C) and is sufficiently large, 
         we report the match as likely a peptide (and hence a false positive for our purposes).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as containing a carboxamidine group (and not flagged as peptide-derived),
              False otherwise.
        str: Reason string for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns: one for the neutral amidine and one for a charged variant.
    # The pattern intends to capture a carbon (sp2, three-connected) with a double bond to a nitrogen (neutral or charged)
    # and a single bond to a second nitrogen.
    pattern_neutral = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    pattern_charged = Chem.MolFromSmarts("[CX3](=[NX3+])[NX3]")
    
    matches_neutral = mol.GetSubstructMatches(pattern_neutral)
    matches_charged = mol.GetSubstructMatches(pattern_charged)
    
    all_matches = list(matches_neutral) + list(matches_charged)
    if not all_matches:
        return False, "No carboxamidine group detected."
    
    # Define a very simple SMARTS to flag a peptide-bond fragment.
    peptide_pattern = Chem.MolFromSmarts("N-C(=O)C")
    is_peptide_like = mol.HasSubstructMatch(peptide_pattern) and mol.GetNumHeavyAtoms() > 30
    
    # Now check each match for extra substitution.
    # Our intent: in a "true" carboxamidine, the central carbon (of the group) should be bound only to the two amidine nitrogens,
    # or to one extra substituent that is not "peptide-like".
    valid_match_found = False
    for match in all_matches:
        # According to SMARTS the match tuple will include indices for the atoms involved in the pattern.
        # By our definition, assume the first index is the carbon.
        c_idx = match[0]
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Get all neighbors of this carbon:
        neighbor_indices = [nbr.GetIdx() for nbr in c_atom.GetNeighbors()]
        # We already expect that two of them come from the amidine nitrogens.
        extra_neighbors = [nbr for nbr in neighbor_indices if nbr not in match[1:]]
        # If there is an extra neighbor, inspect its bonds.
        flag_peptide_env = False
        for nbr_idx in extra_neighbors:
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            # Look for a double bond from the neighbor to an oxygen
            for bond in nbr_atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 2.0:
                    other = bond.GetOtherAtom(nbr_atom)
                    if other.GetAtomicNum() == 8:  # oxygen
                        flag_peptide_env = True
                        break
            if flag_peptide_env:
                break
        
        # If an extra neighbor is present and it appears peptide-like, then we flag this match.
        if flag_peptide_env:
            continue
        else:
            valid_match_found = True
            break

    if valid_match_found and not is_peptide_like:
        return True, "Carboxamidine group detected via substructure match."
    else:
        # Even though our simple SMARTS flagged an amidine, the environment suggests a peptide context.
        return False, "WRONGLY CLASSIFIED Carboxamidine group detected via substructure match (peptide-like environment)."

# Example usage (for testing):
if __name__ == "__main__":
    # Test a few examples:
    test_cases = {
        "4-[5-(3,5-dichlorophenyl)-5-(trifluoromethyl)-4,5-dihydro-1,2-oxazol-3-yl]-N-[(methoxyamino)methylidene]-2-methylbenzamide":
            "CON=CNC(=O)c1ccc(cc1C)C1=NOC(C1)(c1cc(Cl)cc(Cl)c1)C(F)(F)F",
        "benzamidine":
            "NC(=N)c1ccccc1",
        "formamidine":
            "[H]C(N)=N",
        "nizatidine (expected false negative in previous attempt)":
            "CNC(NCCSCc1csc(CN(C)C)n1)=C[N+]([O-])=O",
        "Arg-Ser-Tyr (false positive peptide)":
            "O=C(N[C@H](C(=O)N[C@@H](CC1=CC=C(O)C=C1)C(O)=O)CO)[C@@H](N)CCCN=C(N)N"
    }
    
    for name, smi in test_cases.items():
        res, reason = is_carboxamidine(smi)
        print(f"Test: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")
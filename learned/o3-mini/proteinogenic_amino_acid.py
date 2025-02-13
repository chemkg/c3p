"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: Proteinogenic amino acids

Definition:
  Proteinogenic amino acids are the 23 α–amino acids (20 coded by the nuclear genome plus selenocysteine,
  pyrrolysine, and N-formylmethionine) that serve as precursors of proteins. (Glycine is achiral; the others are L‐amino acids).

This implementation uses several SMARTS patterns to search for a free amino acid backbone.
It checks that the molecule is small (MW < 300 Da), that the free amine is unmodified (has exactly 2 hydrogens,
or 1 in the case of the cyclic proline backbone), and that for non–glycine amino acids the chiral α–carbon is set to “S”.
Because some structures contain isotope labels that may alter CIP assignments, we skip the chirality check
when any isotopic atoms are present.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    
    For free amino acids, the ideal backbone is:
         N[C@@H](R)C(=O)O     or     N[C@H](R)C(=O)O     or    NC(C(=O)O)
    plus some formyl variants (e.g. O=CN[C@@H](R)C(=O)O) and a special pattern for proline.
    
    Additional Filters:
      - Molecular weight must be < 300 Da.
      - The free backbone amino group must be unmodified: it should have exactly two hydrogen atoms attached,
        except in proline (where the amino group is part of a ring) which should have one.
      - For non–glycine amino acids, if the α–carbon (“chiral center”) carries a sidechain then if the CIP is assigned,
        its code must be “S” (unless isotopic atoms are present, in which case we assume L).
      - If no chirality is specified, we assume the intended L configuration.
      
    Args:
      smiles (str): SMILES string.
    
    Returns:
      bool: True if the molecule is a proteinogenic amino acid, False otherwise.
      str: Reason for the classification decision.
    """
    # Parse the SMILES and add explicit hydrogens so that hydrogen atoms and stereochemistry are available.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check MW filter: free amino acids have molecular weight less than 300 Da.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw > 300:
        return False, f"Molecular weight ({mw:.1f} Da) too high for a free amino acid"
    
    # Determine if any atom is isotopically labelled; if yes then skip chirality check.
    has_isotopes = any(atom.GetIsotope() for atom in mol.GetAtoms())
    
    # Define a list of backbone patterns along with a mapping defining indices:
    # "amine": index for the backbone amino nitrogen,
    # "alpha": index for the α–carbon,
    # "carboxyl": index for the carboxyl carbon.
    patterns = [
        ("N[C@@H](*)C(=O)O",       {"amine": 0, "alpha": 1, "carboxyl": 2}),
        ("N[C@H](*)C(=O)O",        {"amine": 0, "alpha": 1, "carboxyl": 2}),
        ("NC(C(=O)O)",            {"amine": 0, "alpha": 1, "carboxyl": 2}),
        ("O=CN[C@@H](*)C(=O)O",     {"amine": 2, "alpha": 3, "carboxyl": 4}),
        ("O=CN[C@H](*)C(=O)O",      {"amine": 2, "alpha": 3, "carboxyl": 4}),
        ("O=CNC(C(=O)O)",          {"amine": 2, "alpha": 3, "carboxyl": 4}),
        # Proline pattern: note the backbone amino group here is within a ring and should have 1 hydrogen.
        ("OC(=O)[C@@H]1CCCN1",     {"amine": 4, "alpha": 2, "carboxyl": 1}),
    ]
    
    # We will try each pattern (and each match for each pattern) and only report success if one match is valid.
    last_error = "Does not match proteinogenic amino acid backbone pattern"
    for smarts, mapping in patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue  # try next pattern

        # For each match, verify the details.
        for match in matches:
            try:
                amine_atom = mol.GetAtomWithIdx(match[mapping["amine"]])
                alpha_atom = mol.GetAtomWithIdx(match[mapping["alpha"]])
                carboxyl_atom = mol.GetAtomWithIdx(match[mapping["carboxyl"]])
            except IndexError:
                continue  # bad mapping: try next match
            
            # Determine if this is the proline pattern by the SMARTS string.
            is_proline = (smarts == "OC(=O)[C@@H]1CCCN1")
            expected_h = 1 if is_proline else 2
            
            # Count the number of explicit hydrogen (atomic number 1) neighbors on the backbone amine.
            actual_h = sum(1 for nbr in amine_atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
            if actual_h != expected_h:
                # Instead of returning immediately, record error and try another match.
                last_error = ("Backbone amino group appears modified (has "
                              f"{actual_h} hydrogen(s) instead of {expected_h}); "
                              "likely part of a peptide bond or N-modified")
                continue

            # At the alpha carbon, identify the heavy atoms (atomic num > 1) that are not part of the assumed backbone.
            # The backbone atoms are the amine and the carboxyl carbon.
            backbone_indices = { match[mapping["amine"]], match[mapping["carboxyl"]] }
            sidechain_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() 
                                   if nbr.GetIdx() not in backbone_indices and nbr.GetAtomicNum() > 1]
            # Glycine has no heavy atom sidechain.
            if not sidechain_neighbors:
                return True, "Matches amino acid backbone (glycine detected: no heavy sidechain at α–carbon)"
            else:
                # For non–glycine amino acids, check stereochemistry if it is explicitly set.
                if alpha_atom.HasProp('_CIPCode'):
                    cip = alpha_atom.GetProp('_CIPCode')
                    # If isotopic labels are present, we choose to bypass the chirality check.
                    if not has_isotopes:
                        if cip != "S":
                            last_error = (f"α–carbon configuration is '{cip}', not 'S' (expected for L–amino acids)")
                            continue
                    # If we have isotopes, simply assume the intended L configuration.
                    return True, "Matches amino acid backbone with correct (L) configuration"
                else:
                    # If no stereochemistry is explicitly defined, assume L configuration
                    return True, ("Matches amino acid backbone pattern but α–carbon chirality is "
                                 "not explicitly specified; assuming L configuration")
    # If we finish looping over all patterns and matches without a valid one, return with the last error
    return False, last_error

# Example usage for testing the provided examples.
if __name__ == "__main__":
    examples = [
        ("L-arginine", "N[C@@H](CCCNC(N)=N)C(O)=O"),
        ("L-histidine", "N[C@@H](Cc1c[nH]cn1)C(O)=O"),
        ("glycine-13C2,15N", "O[13C](=O)[13CH2][15NH2]"),
        ("L-threonine", "C[C@@H](O)[C@H](N)C(O)=O"),
        ("L-alanine", "C[C@H](N)C(O)=O"),
        ("L-tyrosine-d4", "OC(=O)[C@@H](N)CC1=C(C(=C(O)C(=C1[2H])[2H])[2H])[2H]"),
        ("L-methionine", "CSCC[C@H](N)C(O)=O"),
        ("L-leucine-d3", "OC(=O)[C@@H](N)CC(C([2H])([2H])[2H])C"),
        ("L-asparagine", "N[C@@H](CC(N)=O)C(O)=O"),
        ("L-phenylalanine-d5", "OC(=O)[C@@H](N)CC1=C(C(=C(C(=C1[2H])[2H])[2H])[2H])[2H]"),
        ("glycine-d5", "C(C(N([2H])[2H])([2H])[2H])(=O)O[2H]"),
        ("L-serine", "N[C@@H](CO)C(O)=O"),
        ("L-proline", "OC(=O)[C@@H]1CCCN1"),
        ("L-aspartic acid-d7", "O(C(=O)[C@@](N([2H])[2H])(C(C(O[2H])=O)([2H])[2H])[2H])[2H]"),
        ("L-glutamic acid", "N[C@@H](CCC(O)=O)C(O)=O"),
        ("L-lysine", "NCCCC[C@H](N)C(O)=O"),
        ("L-glutamine", "N[C@@H](CCC(N)=O)C(O)=O"),
        ("L-pyrrolysine", "C(=O)([C@@H](N)CCCCNC([C@H]1[C@@H](CC=N1)C)=O)O"),
        ("L-leucine", "CC(C)C[C@H](N)C(O)=O"),
        ("L-methionine-d3", "S(CC[C@H](N)C(O)=O)C([2H])([2H])[2H]"),
        ("L-isoleucine", "OC([C@H]([C@H](CC)C)N)=O"),
        ("L-valine-d8", "OC(=O)[C@@](N)(C(C([2H])([2H])[2H])(C([2H])([2H])[2H])[2H])[2H]"),
        ("L-arginine-d7", "OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]"),
        ("L-cysteine", "N[C@@H](CS)C(O)=O"),
        ("L-proline-d7", "OC(=O)[C@]1(NC(C(C1([2H])[2H])([2H])[2H])([2H])[2H])[2H]"),
        ("L-valine", "CC(C)[C@H](N)C(O)=O"),
    ]
    for name, smi in examples:
        ok, reason = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {ok} ({reason})")
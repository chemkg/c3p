"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: Proteinogenic amino acids

Definition:
  Proteinogenic amino acids are the 23 α–amino acids (20 coded by the nuclear genome plus selenocysteine,
  pyrrolysine, and N-formylmethionine) that serve as protein precursors. Apart from glycine (which is achiral),
  they have an L–configuration. 

This implementation uses SMARTS patterns to search for the free amino acid backbone.
In addition to matching the connectivity (an amino group, a chiral α–carbon, and a carboxylic acid or a formylated variant),
we perform two extra filters:
  • The molecule’s overall molecular weight must be consistent with a single amino acid (< 300 Da).
  • The “free” amino group is checked to ensure that it is not modified into an amide (the backbone N should carry 
    2 hydrogens for a primary amine – except in proline, where it is part of a cycle and carries 1 hydrogen).

For non–glycine amino acids (i.e. those with a sidechain), if stereochemistry (CIP) is available it must be “S”.
If no CIP label is set, we allow the molecule but issue a relaxed result (assuming the L configuration).
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    
    For free amino acids the backbone is typically:
       N[C@@H](R)C(=O)O    or     N[C@H](R)C(=O)O    or      NC(C(=O)O)
    plus formyl variants (e.g., O=CN[C@@H](R)C(=O)O).
    For proline the backbone pattern (with the amino group in the ring) is
       OC(=O)[C@@H]1CCCN1

    In addition:
      • For non–glycine amino acids the free amino group must be unmodified: it should have two hydrogens (or one for proline).
      • The α–carbon: if it carries any non–hydrogen sidechain then (a) if it is chiral it must have CIP “S”; 
        if chirality is not explicitly present we assume L configuration.
      • The molecule must be small (MW < 300 Da) to avoid peptides.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a (free) proteinogenic amino acid, False otherwise.
        str: Explanation/reason for the classification decision.
    """
    # Parse SMILES and add explicit hydrogens so that stereochemistry and H counts are available.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Filter out molecules too heavy to be single amino acids.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw > 300:
        return False, f"Molecular weight ({mw:.1f} Da) too high for a free amino acid"
    
    # Define SMARTS patterns with mapping for:
    #    "alpha": index of α–carbon,
    #    "amine": index of the backbone nitrogen,
    #    "carboxyl": index of the carboxyl carbon.
    # Note: The first six patterns cover the typical backbone (free amine),
    # while the last one is for proline.
    patterns = [
        ("N[C@@H](*)C(=O)O", {"alpha": 1, "amine": 0, "carboxyl": 2}),
        ("N[C@H](*)C(=O)O", {"alpha": 1, "amine": 0, "carboxyl": 2}),
        ("NC(C(=O)O)", {"alpha": 1, "amine": 0, "carboxyl": 2}),
        ("O=CN[C@@H](*)C(=O)O", {"alpha": 3, "amine": 2, "carboxyl": 4}),
        ("O=CN[C@H](*)C(=O)O", {"alpha": 3, "amine": 2, "carboxyl": 4}),
        ("O=CNC(C(=O)O)", {"alpha": 3, "amine": 2, "carboxyl": 4}),
        # Proline: note that in proline the amino group is in a ring and typically has 1 hydrogen.
        ("OC(=O)[C@@H]1CCCN1", {"alpha": 2, "amine": 6, "carboxyl": 1}),
    ]
    
    # Check each pattern for a substructure match.
    for smarts, mapping in patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # skip if pattern is invalid (should not occur)
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue  # try next pattern if no match
        
        # For each match (if multiple, we check until one appears that passes additional filters)
        for match in matches:
            # Make sure the match tuple is long enough.
            try:
                alpha_atom = mol.GetAtomWithIdx(match[mapping["alpha"]])
                amine_atom = mol.GetAtomWithIdx(match[mapping["amine"]])
                carboxyl_atom = mol.GetAtomWithIdx(match[mapping["carboxyl"]])
            except IndexError:
                continue
            
            # Check that the free amino (backbone) N is not amidated.
            # For a free primary amine, after adding hydrogens,
            # the number of hydrogens should be 2. (For proline, it should be 1.)
            expected_h = 1 if "1" in smarts and "CCCN1" in smarts else 2
            actual_h = amine_atom.GetTotalNumHs()
            if actual_h != expected_h:
                return False, ("Backbone amino group appears modified (has "
                               f"{actual_h} hydrogen(s) instead of {expected_h}); "
                               "likely part of a peptide bond")
            
            # Determine if the α–carbon is glycine (no heavy atom on the sidechain).
            # Exclude neighbors that are already part of the backbone.
            backbone_idxs = { match[mapping["amine"]], match[mapping["carboxyl"]] }
            sidechain_neighbors = [nbr for nbr in alpha_atom.GetNeighbors()
                                   if nbr.GetIdx() not in backbone_idxs and nbr.GetAtomicNum() > 1]
            if not sidechain_neighbors:
                # Glycine: lack of heavy-atom sidechain; chirality is not relevant.
                return True, "Matches amino acid backbone (glycine identified: α–carbon is achiral)"
            else:
                # For non-glycine, if the α–carbon carries any heavy-atom sidechain, it should be chiral.
                if alpha_atom.HasProp('_CIPCode'):
                    cip = alpha_atom.GetProp('_CIPCode')
                    if cip != "S":
                        return False, (f"α–carbon configuration is '{cip}', not 'S' (expected for L–amino acids)")
                    else:
                        return True, "Matches amino acid backbone with correct (L) configuration"
                else:
                    # If stereochemistry is not explicitly set but the connectivity is correct,
                    # we assume that the molecule is intended to be L.
                    return True, ("Matches amino acid backbone pattern but α–carbon chirality is "
                                  "not explicitly specified; assuming L configuration")
    
    # If none of the patterns matched, the molecule does not present a free amino acid backbone.
    return False, "Does not match proteinogenic amino acid backbone pattern"

# Example usage testing some of the provided examples.
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
        ("L-aspartic acid", "N[C@@H](CC(O)=O)C(O)=O"),
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
        ok, msg = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {ok} ({msg})")
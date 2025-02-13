"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: Proteinogenic amino acids

Definition:
  Proteinogenic amino acids are the 23 α–amino acids (20 coded by the nuclear genome plus selenocysteine,
  pyrrolysine, and N-formylmethionine) that serve as precursors of proteins. (Glycine is achiral; the others are L‐amino acids).

This implementation uses SMARTS patterns to search for a free amino acid backbone.
In addition to matching the connectivity (a free amino group, a chiral α–carbon, and a carboxylic acid or formyl variant),
we apply two extra filters:
  • The molecule’s molecular weight must be consistent with a single amino acid (< 300 Da).
  • The backbone amino nitrogen must have the expected number of hydrogen atoms attached (2 for free primary amines,
    or 1 for proline where the amino group is part of a ring).

Note: Instead of using GetTotalNumHs, which may return 0 in some cases, we count explicit hydrogen neighbors.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    
    For free amino acids, the backbone is typically:
         N[C@@H](R)C(=O)O    or    N[C@H](R)C(=O)O    or     NC(C(=O)O)
    plus some formyl variants (e.g., O=CN[C@@H](R)C(=O)O).
    For proline the pattern appears as:
         OC(=O)[C@@H]1CCCN1

    Additional Filters:
      - The molecule must be small (MW < 300 Da).
      - The “free” backbone amino group must be unmodified (i.e. it should have two hydrogen atoms attached,
        except in proline, which should have one).
      - For non–glycine amino acids, if the α–carbon carries a sidechain then (if chiral) it must have configuration “S”.
        If no stereochemistry is explicitly set we assume it is L.
    
    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if the molecule is classified as a proteinogenic amino acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES and add explicit hydrogens so that hydrogen atoms and stereochemistry are available.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Filter out molecules that are too heavy to be a single amino acid.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw > 300:
        return False, f"Molecular weight ({mw:.1f} Da) too high for a free amino acid"
    
    # Define SMARTS patterns and corresponding atom mappings.
    # The mapping indicates the indices for:
    #   "amine": backbone amino nitrogen,
    #   "alpha": alpha carbon,
    #   "carboxyl": carboxyl carbon.
    # The first several patterns are for typical amino acid backbones;
    # the last one is for proline.
    patterns = [
        ("N[C@@H](*)C(=O)O", {"amine": 0, "alpha": 1, "carboxyl": 2}),
        ("N[C@H](*)C(=O)O",  {"amine": 0, "alpha": 1, "carboxyl": 2}),
        ("NC(C(=O)O)",       {"amine": 0, "alpha": 1, "carboxyl": 2}),
        ("O=CN[C@@H](*)C(=O)O", {"amine": 2, "alpha": 3, "carboxyl": 4}),
        ("O=CN[C@H](*)C(=O)O",  {"amine": 2, "alpha": 3, "carboxyl": 4}),
        ("O=CNC(C(=O)O)",       {"amine": 2, "alpha": 3, "carboxyl": 4}),
        # Proline pattern – note in proline the amino group is part of the ring so it normally has one hydrogen.
        ("OC(=O)[C@@H]1CCCN1", {"amine": 4, "alpha": 2, "carboxyl": 1}),
    ]
    
    # Loop over patterns to find a matching free amino acid backbone.
    for smarts, mapping in patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # Skip invalid SMARTS pattern (should not occur)
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue  # Try next pattern if no match found
        
        # Evaluate each match in turn.
        for match in matches:
            try:
                amine_atom = mol.GetAtomWithIdx(match[mapping["amine"]])
                alpha_atom = mol.GetAtomWithIdx(match[mapping["alpha"]])
                carboxyl_atom = mol.GetAtomWithIdx(match[mapping["carboxyl"]])
            except IndexError:
                continue  # If match tuple is not as expected, try next match
            
            # Determine expected hydrogen count: proline (our only cycle pattern) should have 1,
            # others should have 2 hydrogens attached to the backbone nitrogen.
            is_proline = smarts == "OC(=O)[C@@H]1CCCN1"
            expected_h = 1 if is_proline else 2
            
            # Instead of GetTotalNumHs(), count the explicit hydrogen atoms
            # bonded to the backbone amine.
            actual_h = sum(1 for nbr in amine_atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
            if actual_h != expected_h:
                return False, ("Backbone amino group appears modified (has "
                               f"{actual_h} hydrogen(s) instead of {expected_h}); "
                               "likely part of a peptide bond")
            
            # Check the sidechain at the alpha carbon.
            # Exclude backbone atoms (amine and carboxyl) to find sidechain heavy atoms.
            backbone_ids = { match[mapping["amine"]], match[mapping["carboxyl"]] }
            sidechain_atoms = [nbr for nbr in alpha_atom.GetNeighbors() 
                               if nbr.GetIdx() not in backbone_ids and nbr.GetAtomicNum() > 1]
            
            # Glycine is the only amino acid with no heavy sidechain.
            if not sidechain_atoms:
                return True, "Matches amino acid backbone (glycine detected: no heavy sidechain at α–carbon)"
            else:
                # For non-glycine amino acids, check chirality if it is explicitly set.
                if alpha_atom.HasProp('_CIPCode'):
                    cip = alpha_atom.GetProp('_CIPCode')
                    if cip != "S":
                        return False, (f"α–carbon configuration is '{cip}', not 'S' (expected for L–amino acids)")
                    else:
                        return True, "Matches amino acid backbone with correct (L) configuration"
                else:
                    # If stereochemistry is not explicitly defined, assume L configuration.
                    return True, ("Matches amino acid backbone pattern but α–carbon chirality is "
                                  "not explicitly specified; assuming L configuration")
    
    # No matching amino acid backbone found
    return False, "Does not match proteinogenic amino acid backbone pattern"

# Example usage for testing some provided examples.
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
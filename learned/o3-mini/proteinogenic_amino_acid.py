"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
#!/usr/bin/env python
"""
Classifies: Proteinogenic amino acids
Definition:
  Any of the 23 alpha‐amino acids that serve as precursors to proteins (the 20 encoded by the nuclear genome of eukaryotes plus selenocysteine, pyrrolysine, 
  and N-formylmethionine). Apart from glycine (which is achiral), all have L configuration.
  
In this implementation we use SMARTS patterns to detect the amino acid backbone.
The patterns include a free amino acid backbone, as well as formylated (e.g. N-formylmethionine) variants.
For a non‐glycine candidate we further check that the chiral α‐carbon has a defined CIP code equal to “S”.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    It uses SMARTS matching to detect a backbone of the form:
      (free amino acid)  N[C@H](*)C(=O)O  or  N[C@@H](*)C(=O)O   or NC(C(=O)O)
    or a formylated variant:
      O=CN[C@H](*)C(=O)O,  O=CN[C@@H](*)C(=O)O or O=CNC(C(=O)O)
    and then checks that the α‐carbon (the carbon attached to both the amine and carboxylic acid)
    is either achiral (glycine) or, if chiral, has CIP code “S”.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a proteinogenic amino acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemical information is computed
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns and their expected index mapping.
    # Each entry is a tuple: (pattern, mapping)
    # mapping is a dict giving the index positions in the match for:
    #   "alpha": the α‐carbon, "amine": the nitrogen bearing group, "carboxyl": the carboxyl carbon.
    patterns = [
        ("N[C@@H](*)C(=O)O", {"alpha": 1, "amine": 0, "carboxyl": 2}),
        ("N[C@H](*)C(=O)O", {"alpha": 1, "amine": 0, "carboxyl": 2}),
        ("NC(C(=O)O)", {"alpha": 1, "amine": 0, "carboxyl": 2}),
        ("O=CN[C@@H](*)C(=O)O", {"alpha": 3, "amine": 2, "carboxyl": 4}),
        ("O=CN[C@H](*)C(=O)O", {"alpha": 3, "amine": 2, "carboxyl": 4}),
        ("O=CNC(C(=O)O)", {"alpha": 3, "amine": 2, "carboxyl": 4}),
    ]
    
    found_match = False
    reason = ""
    for smarts, mapping in patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # skip if SMARTS invalid
        matches = mol.GetSubstructMatches(patt)
        if matches:
            # We take the first match (heuristic)
            match = matches[0]
            alpha_idx = mapping["alpha"]
            amine_idx = mapping["amine"]
            carboxyl_idx = mapping["carboxyl"]
            # Check that the match is long enough
            if alpha_idx >= len(match) or amine_idx >= len(match) or carboxyl_idx >= len(match):
                continue
            alpha_atom = mol.GetAtomWithIdx(match[alpha_idx])
            # Determine if the α‐carbon is glycine (i.e. its sidechain is only hydrogens)
            # Get neighbors of the α‐carbon but exclude the atoms that are part of the backbone (amine and carboxyl)
            backbone_idxs = { match[amine_idx], match[carboxyl_idx] }
            sidechain_nbrs = [nbr for nbr in alpha_atom.GetNeighbors() 
                              if nbr.GetIdx() not in backbone_idxs and nbr.GetAtomicNum() > 1]
            if not sidechain_nbrs:
                # No heavy-atom sidechain => glycine. Glycine is allowed.
                found_match = True
                reason = "Matches amino acid backbone (glycine identified: α‐carbon is achiral)."
                break
            else:
                # For non-glycine amino acids, the α‐carbon should be chiral with defined stereochemistry.
                if not alpha_atom.HasProp('_ChiralityPossible'):
                    return False, "α‐carbon is expected to be chiral but is not marked as such"
                # Check if CIP code is computed; if not, try to assign it.
                try:
                    cip = alpha_atom.GetProp('_CIPCode')
                except KeyError:
                    return False, "α‐carbon chiral center lacks defined configuration"
                # For L-amino acids (except cysteine complications) we expect CIP code "S"
                if cip != "S":
                    return False, f"α‐carbon configuration is '{cip}', not 'S' (expected for L-amino acids)"
                found_match = True
                reason = "Matches free or formylated amino acid backbone pattern with correct (L) configuration"
                break

    if found_match:
        return True, reason
    else:
        return False, "Does not match proteinogenic amino acid backbone pattern"

# Example usage (you can remove or comment these lines when deploying as a module)
if __name__ == "__main__":
    # Test examples from the prompt:
    examples = [
        ("L-arginine", "N[C@@H](CCCNC(N)=N)C(O)=O"),
        ("L-histidine", "N[C@@H](Cc1c[nH]cn1)C(O)=O"),
        ("glycine-13C2,15N", "O[13C](=O)[13CH2][15NH2]"),
        ("L-threonine", "C[C@@H](O)[C@H](N)C(O)=O"),
        ("L-alanine", "C[C@H](N)C(O)=O"),
        ("L-tyrosine-d4", "OC(=O)[C@@H](N)CC1=C(C(=C(O)C(=C1[2H])[2H])[2H])[2H"),
        ("L-methionine", "CSCC[C@H](N)C(O)=O"),
    ]
    for name, smi in examples:
        ok, msg = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {ok} ({msg})")
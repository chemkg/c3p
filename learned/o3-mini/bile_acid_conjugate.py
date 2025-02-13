"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugate (any bile acid conjugated to a group that increases hydrophilicity/charge)
Based on the definition:
  “Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule.
   Molecules used for conjugation are: glycine, taurine (and other amino acids); sulfuric acid (sulfate);
   glucuronic acid (glucuronate); glucose and other uncharged sugars; and coenzyme A.”
Note: This implementation uses heuristic SMARTS patterns to first detect a steroid core (bile acid backbone)
and then search for common conjugation moieties.
"""
from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a bile acid conjugate.
    The algorithm first looks for the typical steroid (cholanic) nucleus and then for
    a conjugation functional group (such as an amide bond as in glycine/other amino acid conjugates,
    taurine-like groups, sulfate, glucuronate, or a sugar moiety).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Check for a bile acid (steroid) core ---
    # Here we use a simplified SMARTS pattern to represent a cholane-like core.
    # Note: This pattern is heuristic and might not catch every bile acid variant.
    steroid_smarts = "[C@]1(C)CC[C@H]2[C@@H]3CC[C@]4(C)[C@@H](CC[C@]43)[C@H]12"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus (cholanic acid core) detected"
    
    # --- Step 2. Check for conjugation groups ---
    # We define a set of SMARTS representing groups used for conjugation.
    # For example, an amide linkage (as seen with glycine or other amino acids),
    # taurine conjugation (with a characteristic NCCS(=O)(=O)[O-] motif or neutral variant),
    # sulfate groups, a glucuronate motif or a sugar moiety.
    conjugation_smarts = {
        "amide bond (e.g. glycine conjugation)": "C(=O)N",  # very generic amide bond
        "taurine conjugation (anionic)": "NCCS(=O)(=O)[O-]",
        "taurine conjugation (neutral)": "NCCS(=O)(=O)O",
        "sulfate conjugation": "S(=O)(=O)[O-]",
        "glucuronate": "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)O1",
        "sugar (e.g. uncharged glucose)": "OC1OC(O)C(O)C(O)C1O"
        # Additional patterns (e.g. for coenzyme A) could be added as needed.
    }
    
    # Check if at least one conjugation pattern matches.
    conjugation_found = False
    matched_conjugations = []
    for desc, smarts in conjugation_smarts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            conjugation_found = True
            matched_conjugations.append(desc)
    
    if not conjugation_found:
        return False, "No recognized conjugation group detected"
    
    # If both a steroid core and at least one conjugation group are detected, classify as a bile acid conjugate.
    return True, ("Steroid nucleus detected with conjugation groups: " +
                  ", ".join(matched_conjugations))
    
# Example usage (uncomment to test):
# smiles_example = "S(O)(=O)(=O)CCNC(=O)[C@@H](CCC[C@H](C1[C@@]2(C(C3C([C@@]4([C@](CC3O)(C[C@H](O)CC4)[H])C)CC2)CC1)C)C)C"
# result, reason = is_bile_acid_conjugate(smiles_example)
# print(result, reason)
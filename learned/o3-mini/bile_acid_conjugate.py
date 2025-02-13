"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugate – any bile acid (a molecule with a cholanic/steroid core)
conjugated to a functional group that increases hydrophilicity or charge.
Based on the definition:
  “Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule.
   Molecules used for conjugation are: glycine, taurine (and other amino acids); sulfuric acid (sulfate);
   glucuronic acid (glucuronate); glucose and other uncharged sugars; and coenzyme A.”
This implementation uses heuristic SMARTS patterns to detect a steroid (cholanic) nucleus and then
to detect common conjugation moieties.
"""
from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a bile acid conjugate.
    The algorithm uses a relaxed SMARTS pattern to identify a cholanic (steroid) nucleus,
    and then looks for at least one conjugation group from a list of common patterns.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Check for a bile acid (steroid) core ---
    # Previous attempt used a strict, chiral SMARTS pattern that missed many bile acids.
    # Here we use a relaxed (achiral) SMARTS pattern for a fused tetracyclic system, 
    # which is common in bile acid (cholanic) cores.
    steroid_smarts = "C1CCC2C3CCC4C1C2C34"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Error creating steroid SMARTS pattern"
    
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus (cholanic acid core) detected"
    
    # --- Step 2. Check for conjugation groups ---
    # Define several SMARTS patterns for conjugation groups.
    # These include a generic amide (which may capture glycine and other amino acid conjugates),
    # taurine conjugation (as either anionic or neutral form),
    # sulfate groups, glucuronate pattern, sugar moieties, and others.
    conjugation_smarts = {
        "amide bond (e.g. amino acid conjugation)": "C(=O)N",  # generic amide bond
        "taurine conjugation (anionic)": "NCCS(=O)(=O)[O-]",
        "taurine conjugation (neutral)": "NCCS(=O)(=O)O",
        "sulfate conjugation": "S(=O)(=O)[O-]",
        "glucuronate": "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)O1",
        "sugar (e.g. uncharged glucose)": "OC1OC(O)C(O)C(O)C1O"
        # Additional patterns (e.g., coenzyme A) can be added here.
    }
    
    found_conjugation = False
    matched_groups = []
    for desc, smarts in conjugation_smarts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_conjugation = True
            matched_groups.append(desc)
    if not found_conjugation:
        return False, "No recognized conjugation group detected"
    
    # If both a steroid nucleus and at least one conjugation group are found, classify as bile acid conjugate.
    return True, "Steroid nucleus detected with conjugation group(s): " + ", ".join(matched_groups)

# Example usage:
if __name__ == "__main__":
    # List of example SMILES strings (e.g., taurocholic acid, taurochenodeoxycholic acid, etc.)
    smiles_examples = [
        "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCC(=O)NCCS(O)(=O)=O",  # taurocholic acid
        "[C@]12([C@]([C@]3([C@]([C@@H]([C@@H](CCC(NCCS(=O)(=O)O)=O)C)CC3)(C)CC1)[H])([C@H]([C@H](O)[C@]4([C@@]2(CC[C@H](C4)O)C)[H])O)[H])[H]"  # tauro-omega-muricholic acid
    ]
    for s in smiles_examples:
        result, reason = is_bile_acid_conjugate(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("----------")
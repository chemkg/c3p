"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugate
Definition:
  Any bile acid (i.e. a molecule with a cholanic/steroid core – typically a fused system
  of three six-membered rings and one five-membered ring) conjugated to a functional group
  that gives additional hydrophilicity or charge (e.g. glycine, taurine, sulfate, glucuronate,
  sugars, or coenzyme A).
  
This implementation uses a heuristic method:
  1. Detect the steroid core by examining the ring systems.
  2. Look for one or more conjugation groups using common SMARTS patterns.
"""

from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a bile acid conjugate.
    It uses a heuristic for the steroid (cholanic) nucleus based on ring counts 
    (expecting 3 six-membered rings and 1 five-membered ring) and then searches 
    for common conjugation groups using SMARTS patterns.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a bile acid conjugate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1. Check for a steroid (cholanic acid) nucleus ---
    # Bile acids have a tetracyclic fused ring system: three cyclohexane rings and one cyclopentane ring.
    # We use RDKit's ring information as a heuristic.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    if len(rings) < 4:
        return False, "Less than 4 rings detected – steroid nucleus likely missing"
    
    # Count how many rings have 5 atoms and how many have 6 atoms.
    count_5 = sum(1 for ring in rings if len(ring) == 5)
    count_6 = sum(1 for ring in rings if len(ring) == 6)
    if count_5 < 1 or count_6 < 3:
        return False, ("Detected rings do not match the expected pattern for a steroid nucleus "
                       "(expect at least one 5-membered ring and three 6-membered rings)")
    
    # --- Step 2. Check for conjugation groups ---
    # Define a list of SMARTS patterns for functional groups used for bile acid conjugation.
    # These include generic amide bonds (common for amino acid conjugates), taurine and sulfate conjugates,
    # glucuronate, and sugar moieties.
    conjugation_smarts = {
        "generic amide (e.g. amino acid conjugation)": "C(=O)N",
        "taurine conjugation (anionic)": "NCCS(=O)(=O)[O-]",
        "taurine conjugation (neutral)": "NCCS(=O)(=O)O",
        "sulfate conjugation": "S(=O)(=O)[O-]",
        "glucuronate": "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-]1",
        "sugar (e.g. glucose)": "OC1OC(O)C(O)C(O)C1O"
        # Note: Additional patterns (e.g., for coenzyme A) can be included here as needed.
    }
    
    found_conjugation = False
    matched_groups = []
    for desc, smart in conjugation_smarts.items():
        pattern = Chem.MolFromSmarts(smart)
        if pattern is None:
            continue  # Skip invalid SMARTS
        if mol.HasSubstructMatch(pattern):
            found_conjugation = True
            matched_groups.append(desc)
            
    if not found_conjugation:
        return False, "No recognized conjugation group detected"
    
    # If both a steroid nucleus and at least one conjugation group are found,
    # we classify the molecule as a bile acid conjugate.
    reason = ("Steroid nucleus detected (heuristic: >=1 five-membered and >=3 six-membered rings) "
              "with conjugation group(s): " + ", ".join(matched_groups))
    return True, reason

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example SMILES from the provided list (taurocholic acid is given as one example)
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
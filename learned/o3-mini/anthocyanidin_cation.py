"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin Cation (flavylium‐type aglycon)
Definition: Any organic cation that is an aglycon of anthocyanin cation;
they are oxygenated derivatives of flavylium (2-phenylchromenylium).

This version improves on the previous code by:
  - Desalting the molecule (breaking into fragments).
  - Checking that the fragment’s net formal charge is positive.
  - Matching a flavylium-type core using a SMARTS pattern.
  - Rejecting fragments that show evidence of glycosylation (sugar rings)
    or have too many rings (the aglycon should have only the fused system).
  
Heuristic criteria:
  - A fragment must have a net positive formal charge.
  - It must contain a flavylium core defined by the SMARTS pattern:
       c1ccc2c(c1)[o+][c]([c]3ccccc3)cc2
    which looks for a six-membered ring containing an oxonium,
    fused to a benzene ring and bearing an attached phenyl substituent.
  - It must not contain any evident sugar moieties (using a simple sugar SMARTS),
    and the number of rings in the fragment should be 3 (the core rings).
    
Note: This heuristic is not perfect but improves on the previous version by filtering
out many glycosylated anthocyanins.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines whether a molecule (or one of its fragments) qualifies as an anthocyanidin cation.
    
    Steps:
      1. Parse and sanitize the input molecule.
      2. Break it into fragments so that counterions do not hide the positive ion.
      3. For each fragment:
           - Check that the net formal charge is positive.
           - Remove hydrogens for substructure matching.
           - Check for a flavylium core using a SMARTS pattern.
           - Reject fragments that appear glycosylated (i.e. have a sugar substructure)
             or that have more than 3 rings.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if at least one fragment is classified as an anthocyanidin cation.
      str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        # Add hydrogens and sanitize to ensure proper aromaticity.
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # Break molecule into fragments (to remove salts/counterions).
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    # Define a SMARTS pattern for a flavylium-type core.
    # The pattern looks for a bicyclic core:
    #  - a benzene ring fused to a six-membered oxygen-containing ring;
    #  - the oxygen in that ring carries a formal +1 charge; and
    #  - one of the ring carbons carries a benzene substituent.
    flavylium_smarts = "c1ccc2c(c1)[o+][c]([c]3ccccc3)cc2"
    flavylium_pattern = Chem.MolFromSmarts(flavylium_smarts)
    if flavylium_pattern is None:
        return False, "Flavylium SMARTS pattern failed to compile"
    
    # Define a simple SMARTS pattern to catch common sugar (pyranose) rings.
    # This pattern looks for a 6-membered oxygen-containing ring with multiple hydroxyls.
    sugar_smarts = "[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    
    # Check each fragment for the flavylium-type core.
    for frag in frags:
        # Check that the fragment has positive net charge.
        if frag.GetFormalCharge() <= 0:
            continue  # Skip if not cationic
        
        # Optionally remove hydrogens for substructure matching.
        frag_noH = Chem.RemoveHs(frag)
        
        # Require that the flavylium core SMARTS is matched.
        if not frag_noH.HasSubstructMatch(flavylium_pattern):
            continue

        # Check ring count: the aglycon should possess the expected 3 rings.
        ri = frag_noH.GetRingInfo()
        if ri.NumRings() > 3:
            # Likely extra rings from glycosides – skip as not an aglycon.
            continue

        # If a sugar substructure is found, then this fragment is glycosylated.
        if sugar_pattern is not None and frag_noH.HasSubstructMatch(sugar_pattern):
            continue

        # (Optional) enforce an upper bound on molecular weight (the aglycon is usually <500 Da).
        mw = AllChem.CalcExactMolWt(frag_noH)
        if mw > 500:
            continue

        # If all criteria pass, return True.
        return True, ("Fragment matches flavylium-type (2-phenylchromenylium) core with positive net charge, "
                      "has 3 rings and no sugar substructure detected – likely an anthocyanidin cation.")
    
    return False, "No fragment met the criteria for an anthocyanidin cation (flavylium aglycon with no glycosides)"

# Example usage:
if __name__ == '__main__':
    # Testing with one of the provided examples (e.g., hirsutidin, which should be classified as True)
    test_smiles = "COc1cc(O)c2cc(O)c([o+]c2c1)-c1cc(OC)c(O)c(OC)c1"
    result, reason = is_anthocyanidin_cation(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
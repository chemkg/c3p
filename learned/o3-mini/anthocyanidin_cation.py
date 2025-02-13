"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin Cation (flavylium‐type aglycon)
Definition: Any organic cation that is an aglycon of anthocyanin cation;
they are oxygenated derivatives of flavylium (2-phenylchromenylium).

This program improves on the previous version by:
  - Calculating the net formal charge manually (summing over atoms), since 
    mol objects do not have a GetFormalCharge method.
  - Desalting the molecule (breaking into fragments).
  - Checking that the fragment’s net formal charge is positive.
  - Matching a flavylium core using a SMARTS pattern.
  - Rejecting fragments that show evidence of glycosylation (sugar rings)
    or have extra rings beyond the expected three fused rings in the aglycon.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines whether a molecule (or one of its fragments) qualifies as an anthocyanidin cation.
    
    Heuristic steps:
      1. Parse and sanitize the input molecule.
      2. Break molecule into fragments to remove salts/counterions.
      3. For each fragment:
           - Calculate the net formal charge manually by summing atom charges.
           - Remove hydrogens for substructure matching.
           - Check for a flavylium core using a SMARTS pattern.
           - Reject fragments that appear glycosylated (using a sugar SMARTS)
             or that have too many rings (more than 3).
           - Optionally, filter by molecular weight (aglycon is usually <500 Da).
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if at least one fragment qualifies as an anthocyanidin cation.
      str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        # Add explicit hydrogens and sanitize to ensure correct aromaticity and charge assignment.
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # Break molecule into fragments (to remove salts and counterions).
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    # Define SMARTS pattern for a flavylium-type core.
    # This pattern tries to match a six-membered ring with an oxonium ([o+])
    # that is fused to a benzene and bearing a phenyl substituent.
    flavylium_smarts = "c1ccc2c(c1)[o+][c]([c]3ccccc3)cc2"
    flavylium_pattern = Chem.MolFromSmarts(flavylium_smarts)
    if flavylium_pattern is None:
        return False, "Flavylium SMARTS pattern failed to compile"
    
    # Define a simple SMARTS pattern to detect sugar (pyranose) rings,
    # which commonly appear in glycosylated anthocyanins.
    sugar_smarts = "[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    
    # Analyze each fragment.
    for frag in frags:
        # Compute net formal charge manually by summing the formal charge of each atom.
        net_charge = sum(atom.GetFormalCharge() for atom in frag.GetAtoms())
        if net_charge <= 0:
            continue  # Skip fragments that are not cationic
        
        # Remove hydrogens for substructure matching.
        frag_noH = Chem.RemoveHs(frag)
        
        # Check if the flavylium-type core is present.
        if not frag_noH.HasSubstructMatch(flavylium_pattern):
            continue
        
        # Check that the ring count is appropriate (aglycon should have 3 rings).
        ring_info = frag_noH.GetRingInfo()
        if ring_info.NumRings() > 3:
            # Extra rings likely indicate glycosylation.
            continue
        
        # Reject fragments that contain a sugar substructure.
        if sugar_pattern is not None and frag_noH.HasSubstructMatch(sugar_pattern):
            continue
        
        # Enforce an optional upper bound on molecular weight; aglycon is usually <500 Da.
        mw = AllChem.CalcExactMolWt(frag_noH)
        if mw > 500:
            continue

        # All criteria are met for this fragment.
        return True, ("Fragment has a flavylium-type core with a positive net charge (net_charge=%d), "
                      "contains the expected ring system, and no detected sugar moieties – "
                      "likely an anthocyanidin cation." % net_charge)
    
    return False, "No fragment met the criteria for an anthocyanidin cation (flavylium aglycon with no glycosylation)"

# Example usage:
if __name__ == '__main__':
    # Testing with one of the provided examples (e.g., hirsutidin)
    test_smiles = "COc1cc(O)c2cc(O)c([o+]c2c1)-c1cc(OC)c(O)c(OC)c1"
    result, reason = is_anthocyanidin_cation(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
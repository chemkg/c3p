"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin Cation (flavylium‐type aglycon)
Definition: Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated derivatives of flavylium (2-phenylchromenylium).
Heuristic approach in this version:
  1. Parse the molecule and break it into fragments (to “remove” counterions).
  2. Instead of requiring that the entire molecule is net positively charged, we now check each fragment.
  3. We define a SMARTS pattern for the flavylium core.
     That is, the pattern looks for a six‐membered ring that contains an aromatic oxygen with formal charge +1 (the “oxonium”),
     fused to a benzene ring, and bearing an attached phenyl substituent (for the 2‑phenylchromenylium core).
  4. If at least one fragment shows the core we consider the molecule positive.
Note: This heuristic does not catch every nuance of the chemistry. It may be improved further by adjusting SMARTS.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines whether a molecule (or one of its fragments) qualifies as an anthocyanidin cation.
    
    Heuristic criteria:
      - The molecule (or one fragment) must contain a flavylium-type (2-phenylchromenylium) core,
        which is characterized by a six-membered oxygen-containing (pyran) aromatic ring
        that carries a formal +1 (oxonium) and is fused with at least one aromatic ring.
      - Also, an aromatic (typically phenyl) substituent should be attached on the correct position.
    
    This function first “de-salts” the input (by breaking into fragments) so that counterions do not hide the cationic core.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one fragment is classified as an anthocyanidin cation
        str: Reason for classification
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (Optional) Sanitize, add explicit hydrogens for proper aromaticity assignment.
    try:
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # Break molecule into fragments (to remove counterions that may mask the positive ion)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    # Define a SMARTS pattern for a flavylium-type core.
    # This pattern attempts to capture a bicyclic system:
    #  - A six-membered oxygen-containing aromatic ring (with the oxygen carrying formal charge +1),
    #    fused to a benzene and with an additional aromatic substituent (phenyl) attached to the cationic ring.
    #
    # The SMARTS below is one heuristic choice. It looks for a substructure:
    #   c1ccc2c(c1)[o+][c]([c]3ccccc3)cc2
    #
    # Explanation:
    #   - "c1ccc2c(c1)" : a benzene ring fused to ring "2"
    #   - "c2" continues the second ring, within which
    #   - "[o+]" is our oxygen (formal charge +1) in the ring,
    #   - "[c]([c]3ccccc3)" means that one of the carbons in the ring carries a substituent that is a benzene.
    #   - "cc2" closes the ring.
    #
    # This pattern may be adjusted further.
    flavylium_smarts = "c1ccc2c(c1)[o+][c]([c]3ccccc3)cc2"
    pattern = Chem.MolFromSmarts(flavylium_smarts)
    if pattern is None:
        return False, "Flavylium SMARTS pattern failed to compile"
    
    # Check each fragment for the flavylium-type core.
    for frag in frags:
        # (Optional) Remove hydrogens to match the substructure in case they interfere.
        frag_noH = Chem.RemoveHs(frag)
        if frag_noH.HasSubstructMatch(pattern):
            return True, ("Fragment matches flavylium-type (2-phenylchromenylium) core: "
                          "fused aromatic rings with an oxonium ([o+]) in a six-membered ring "
                          "and an attached phenyl substituent detected.")
    
    return False, "No fragment contained a flavylium (2-phenylchromenylium) core with a positively charged oxygen"

# Example usage:
if __name__ == '__main__':
    # You can test one of the provided examples (e.g., hirsutidin)
    test_smiles = "COc1cc(O)c2cc(O)c([o+]c2c1)-c1cc(OC)c(O)c(OC)c1"
    result, reason = is_anthocyanidin_cation(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
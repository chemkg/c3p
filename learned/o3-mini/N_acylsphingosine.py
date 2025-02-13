"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N-acylsphingosine
Definition: The parent compounds of the ceramide family, composed of sphingosine having an unspecified fatty acyl group attached to the nitrogen.
This version improves on the previous attempt by:
  - Rejecting molecules that contain any ring systems (which are usually part of glycosylations or other modifications).
  - Lowering the minimum total number of carbon atoms to allow smaller acyl groups (e.g. N-acetylsphingosine).
  - Rejecting molecules containing phosphorus (typically phosphorylated derivatives).
  - Checking for both an amide bond (N-acyl linkage) and a sphingosine backbone fragment.
  - Verifying that the atoms matching the sphingosine backbone SMARTS are acyclic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    
    The following criteria are applied:
      1. The SMILES string must be valid.
      2. The overall molecule must be acyclic (no rings). Parent N-acylsphingosines are not glycosylated or cyclic.
      3. The molecule must contain a minimum number of carbon atoms (>=20).
      4. Molecules containing phosphorus atoms are rejected as these are normally phosphorylated derivatives.
      5. The molecule must contain an amide bond (N-acyl linkage).
      6. The molecule must contain a sphingosine backbone fragment. The heuristic SMARTS used is:
             "[CX4;!R]([OX2H])[CX4;!R](CO)[NX3]C(=O)"
         which looks for:
            - a non‐ring sp3 carbon bearing an –OH,
            - followed by a non‐ring sp3 carbon with a CH2OH substituent,
            - followed by an amide nitrogen.
      7. Additionally, the atoms in any match of the sphingosine backbone pattern must all be acyclic.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is determined to be an N-acylsphingosine, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Requirement 1: Reject molecules that contain any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures; likely glycosylated/cyclic derivative"
    
    # Requirement 2: Check total number of carbon atoms (expect at least ~20 for N-acylsphingosine).
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Too few carbons (found {c_count}); typical N-acylsphingosine have at least ~20 carbons"
    
    # Requirement 3: Reject molecules containing phosphorus (P)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely a phosphorylated derivative"
    
    # Requirement 4: Check that the molecule contains an amide bond (N-acyl linkage).
    # A simple SMARTS for an amide bond is used.
    amide_pattern = Chem.MolFromSmarts("N[C](=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (N-acyl linkage) found"
    
    # Requirement 5: Look for a sphingosine backbone fragment.
    # This SMARTS looks for:
    #   - a tetrahedral, non‐ring carbon bearing an –OH,
    #   - a neighboring tetrahedral, non‐ring carbon with a CH2OH,
    #   - and an attached amide nitrogen.
    sphingosine_smarts = "[CX4;!R]([OX2H])[CX4;!R](CO)[NX3]C(=O)"
    sphingosine_pattern = Chem.MolFromSmarts(sphingosine_smarts)
    matches = mol.GetSubstructMatches(sphingosine_pattern)
    if not matches:
        return False, "Sphingosine backbone pattern not found"
    
    # For each match, verify that the atoms in the matched backbone fragment are all acyclic.
    for match in matches:
        if all(not mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            return True, "Contains an acyclic sphingosine backbone with an N-acyl group attached to the nitrogen"
    
    return False, "Sphingosine backbone pattern detected only in ring environments; expected acyclic backbone"

# Example usage
if __name__ == "__main__":
    # Testing with an example: N-(omega-hydroxyoctacosanoyl)sphingosine
    test_smiles = "C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCO)CO"
    result, reason = is_N_acylsphingosine(test_smiles)
    print("N-acylsphingosine:" if result else "Not N-acylsphingosine:", reason)
"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N-acylsphingosine
Definition: The parent compounds of the ceramide family, composed of sphingosine having an unspecified fatty acyl group attached to the nitrogen.
This version improves on the previous attempt by:
  - Removing a global acyclic check and instead checking that the sphingosine backbone fragment (found by a SMARTS pattern) is acyclic.
  - Rejecting molecules that contain phosphorus (typically found in phosphorylated derivatives).
  - Requiring a minimum total number of carbon atoms.
  - Using a SMARTS pattern that heuristically identifies a sphingosine headgroup composed of an sp3 carbon with an -OH, a neighboring sp3 carbon with a CH2OH substituent, and an adjoining amide bond (the N-acyl linkage).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    
    The following criteria are applied:
      1. The SMILES string must be valid.
      2. The overall molecule must contain a sufficiently high number (>22) of carbon atoms.
      3. Molecules containing phosphorus atoms are rejected as these are typically phosphorylated derivatives.
      4. Look for the presence of an amide bond (the N-acyl linkage).
      5. Look for a sphingosine backbone fragment. Here we use a heuristic SMARTS:
             "[CX4;!R]([OX2H])[CX4;!R](CO)[NX3]C(=O)"
         which looks for:
            - a tetrahedral, non‐ring carbon carrying an –OH group,
            - a neighboring tetrahedral, non‐ring carbon with a CH2OH substituent,
            - directly attached to an amide nitrogen (N–C(=O)).
      6. Instead of requiring the entire molecule to be acyclic, we check that the atoms that match the sphingosine backbone pattern are not in any rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is determined to be an N-acylsphingosine, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Requirement 1: Check total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 23:
        return False, f"Too few carbons (found {c_count}); typical N-acylsphingosine have at least ~23 carbons"
    
    # Requirement 2: Reject molecules containing phosphorus (P)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely a phosphorylated derivative"
    
    # Requirement 3: Check that the molecule contains an amide (N-acyl linkage)
    # Use a simple SMARTS for an amide bond: note this is a coarse check.
    amide_pattern = Chem.MolFromSmarts("N[C](=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (N-acyl linkage) found"
    
    # Requirement 4: Look for a sphingosine backbone fragment.
    # This SMARTS requires:
    #   - a non-ring sp3 carbon bearing an -OH ([CX4;!R]([OX2H]))
    #   - followed by a non-ring sp3 carbon with a CH2OH substituent ([CX4;!R](CO))
    #   - followed by a nitrogen (as part of an amide) ([NX3]C(=O))
    sphingosine_smarts = "[CX4;!R]([OX2H])[CX4;!R](CO)[NX3]C(=O)"
    sphingosine_pattern = Chem.MolFromSmarts(sphingosine_smarts)
    matches = mol.GetSubstructMatches(sphingosine_pattern)
    if not matches:
        return False, "Sphingosine backbone pattern not found"
    
    # For each match, verify that the atoms in the matched backbone fragment are all acyclic.
    for match in matches:
        # match is a tuple of atom indices corresponding to the pattern.
        # Check that none of these atoms is marked as being in a ring.
        if all(not mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            return True, "Contains an acyclic sphingosine backbone with an N-acyl group attached to the nitrogen"
    
    return False, "Sphingosine backbone pattern found only in ring environments; expected acyclic backbone"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: N-(omega-hydroxyoctacosanoyl)sphingosine
    test_smiles = "C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCO)CO"
    result, reason = is_N_acylsphingosine(test_smiles)
    print("N-acylsphingosine:" if result else "Not N-acylsphingosine:", reason)
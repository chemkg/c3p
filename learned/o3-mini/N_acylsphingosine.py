"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N-acylsphingosine
Definition: The parent compounds of the ceramide family, composed of sphingosine having an unspecified fatty acyl group attached to the nitrogen.
Improved classification: In addition to verifying an N-acyl (amide) linkage and a sphingosine-like substructure, we now require
that the molecule is acyclic and has at least a threshold number of carbon atoms.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    
    This improved version applies several filters:
      - The SMILES must be valid.
      - The molecule must contain at least one amide bond (N-acyl linkage).
      - The molecule must be acyclic (sphingosine backbones are typically uncyclized).
      - The molecule must have a sufficiently high total number of carbon atoms.
      - The molecule must contain a sphingosine backbone fragment. Here we look for a fragment having:
            * a sp3 carbon bearing an -OH group,
            * a neighboring sp3 carbon with a CH2OH substituent,
            * and the amide bond (N attached to a carbonyl) directly after.
      The heuristic SMARTS used is:
            "[CX4;!R]([OX2H])[CX4;!R](CO)[NX3]C(=O)"
      which should capture many true N-acylsphingosine structures while reducing false positives.
    
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
    
    # Requirement 1: Molecule should be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; expecting an acyclic sphingosine backbone"
    
    # Requirement 2: Check that the molecule has a sufficiently high number of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Too few carbons (found {c_count}); typical N-acylsphingosine have at least ~18 carbons"
    
    # Requirement 3: Check for an amide bond (the N-acyl linkage).
    amide_pattern = Chem.MolFromSmarts("N[C](=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (N-acyl linkage) found"
    
    # Requirement 4: Look for the sphingosine backbone pattern.
    # This heuristic pattern requires:
    #   - an sp3 carbon with an -OH group ([CX4;!R]([OX2H]))
    #   - a neighboring sp3 carbon with a CH2OH substituent ([CX4;!R](CO))
    #   - an amino group attached that leads into a carbonyl (i.e. the N-acyl part) ([NX3]C(=O))
    sphingosine_pattern = Chem.MolFromSmarts("[CX4;!R]([OX2H])[CX4;!R](CO)[NX3]C(=O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Sphingosine backbone pattern not found"
    
    # If all conditions pass, assume it is an N-acylsphingosine.
    return True, "Contains a sphingosine backbone with an N-acyl group attached to the nitrogen"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: N-(omega-hydroxyoctacosanoyl)sphingosine
    test_smiles = "C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCO)CO"
    result, reason = is_N_acylsphingosine(test_smiles)
    print("N-acylsphingosine:" if result else "Not N-acylsphingosine:", reason)
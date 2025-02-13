"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: endocannabinoid
A class of cannabinoids present in mammalian biological fluids and tissues that activate cannabinoid receptors.
This program uses structural heuristics based on common endocannabinoid substructures:
  - Acyl ethanolamides (e.g., anandamide, palmitoyl ethanolamide)
  - Monoacylglycerols (e.g., 2-arachidonoylglycerol)
It also requires a minimum number of carbons to catch the fatty acyl chain feature.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    
    Heuristic criteria:
      1. The molecule must have a minimum number of carbon atoms (here chosen as 16)
         to suggest the presence of a long fatty acid chain.
      2. The molecule must contain at least one of the two key structural motifs:
           a. Acyl ethanolamide motif: represented by the SMARTS "[C](=O)[N]CCO"
              (found in anandamide, palmitoyl ethanolamide, etc.)
           b. Ester linked to a glycerol backbone: identified by the simultaneous presence
              of an ester group ("OC(=O)") and a glycerol signature ("C(CO)CO")
              (found in monoacylglycerol endocannabinoids such as 2-arachidonoylglycerol).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an endocannabinoid according to the heuristics, False otherwise.
        str: Explanation for the decision.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms - endocannabinoids typically have a long carbon chain.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, "Too few carbon atoms to form a long fatty acyl chain characteristic of endocannabinoids"
    
    # Define SMARTS patterns for the motifs.
    # (a) Acyl ethanolamide motif:
    ethanolamide_smarts = "[C](=O)[N]CCO"
    ethanolamide_pattern = Chem.MolFromSmarts(ethanolamide_smarts)
    
    # (b) Monoacylglycerol motif: an ester group plus a glycerol backbone signature.
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    glycerol_smarts = "C(CO)CO"  # a simple pattern representing the glycerol moiety
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    
    # Check for the acyl ethanolamide feature.
    if mol.HasSubstructMatch(ethanolamide_pattern):
        return True, "Contains an acyl ethanolamide moiety, characteristic of many endocannabinoids."
    
    # Check for the monoacylglycerol feature (ester group on a glycerol backbone).
    if mol.HasSubstructMatch(ester_pattern) and mol.HasSubstructMatch(glycerol_pattern):
        return True, "Contains an ester-linked acyl chain on a glycerol backbone, typical of monoacylglycerol endocannabinoids."
    
    return False, "Does not contain characteristic endocannabinoid structural motifs."
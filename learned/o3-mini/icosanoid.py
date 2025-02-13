"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation 
of the three C20 essential fatty acids (EFAs) icosapentaenoic acid (EPA), arachidonic acid (AA) 
and dihomo-gamma-linolenic acid (DGLA).

Heuristic approach:
 - Parse the SMILES using RDKit.
 - Use a regex on the molecular formula (calculated via RDKit) to extract the number of carbons.
   We expect a total carbon count roughly in the range 15-35 (to allow for modifications).
 - Check if the molecule has at least one carboxylic acid group (SMARTS: "[CX3](=O)[OX1H]").
 - Look for one or both of two characteristic substructures:
        • a cyclopentane ring (SMARTS "C1CCCC1"), which is common in many prostaglandins, or 
        • a polyene chain, using the SMARTS pattern "C=C/C=C/C=C/C=C".
If these conditions are met, we classify the molecule as a potential icosanoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import re

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    (Heuristics: expects a molecule derived from a C20 polyunsaturated fatty acid oxidation,
     usually containing a carboxylic acid and either a cyclopentane ring (prostaglandin core)
     or a long polyene chain.)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is predicted to be an icosanoid, False otherwise
        str: Reason for classification decision
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get molecular formula and extract carbon count
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Using regex: look for a pattern like C12 or C (which implies 1)
    m = re.search(r'C(\d+)', formula)
    if m:
        c_count = int(m.group(1))
    else:
        # If no explicit carbon number found, check atoms manually.
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Many icosanoids arise from C20 fatty acids,
    # so we use a lenient window in case of conjugation: typically between 15 and 35 carbons.
    if not (15 <= c_count <= 35):
        return False, f"Carbon count ({c_count}) is not in the expected range (15–35) for an icosanoid"
    
    # Check for carboxylic acid functionality:
    # SMARTS for carboxylic acid: carbonyl (C=O) with OH
    ca_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if not mol.HasSubstructMatch(ca_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for cyclopentane ring (common in prostaglandin cores)
    cyclopentane = Chem.MolFromSmarts("C1CCCC1")
    has_cyclopentane = mol.HasSubstructMatch(cyclopentane)
    
    # Check for a polyene chain: here we demand at least four conjugated double bonds in sequence
    polyene = Chem.MolFromSmarts("C=C/C=C/C=C/C=C")
    has_polyene = mol.HasSubstructMatch(polyene)
    
    # At least one of these motifs is expected.
    if has_cyclopentane:
        core_reason = "Contains a cyclopentane ring characteristic of prostaglandin cores"
    elif has_polyene:
        core_reason = "Contains a polyene chain (multiple conjugated double bonds) typical of oxidized fatty acids"
    else:
        return False, "Neither a cyclopentane ring nor a polyene chain was detected"
    
    # Optionally, one could check the degree of unsaturation
    # (a high number of double bonds is expected in oxidized fatty acid metabolites)
    unsat = rdMolDescriptors.CalcNumAliphaticRings(mol)  # Note: this counts rings, not double bonds
    # Here, we proceed as our main criteria are the substructures.
    
    return True, (f"Carbon count {c_count} is in range, carboxylic acid present, and {core_reason}.")

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example SMILES: 20-hydroxyprostaglandin A1
    smiles_example = "[C@H]1([C@H](C=CC1=O)/C=C/[C@H](CCCCCO)O)CCCCCCC(O)=O"
    result, reason = is_icosanoid(smiles_example)
    print("Is icosanoid:", result)
    print("Reason:", reason)
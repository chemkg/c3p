"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: Long-chain fatty acyl-CoA(4-)
Definition:
  A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups 
  of any long-chain fatty acyl-CoA; major species at pH 7.3.

This classifier checks for:
  - A CoA moiety (via an adenine-like substructure, SMARTS: n1cnc2c(N)ncnc12).
  - A thioester bond connecting the fatty acyl chain to the CoA. 
    Note: In our molecules the bond may appear as either [CX3](=O)[SX1] or its reverse [SX1]C(=O).
  - After "cleaving" the thioester bond, the fatty acyl chain must have at least 12 carbon atoms.
  - The overall molecule must contain at least 4 deprotonated oxygens (O atoms with a -1 formal charge).
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a given molecule (SMILES string) is a long-chain fatty acyl-CoA(4-).

    Criteria:
      1. SMILES must be valid.
      2. Contains a CoA moiety identified by an adenine-like fragment (SMARTS: n1cnc2c(N)ncnc12).
      3. Contains a thioester bond connecting the acyl chain and CoA. 
         We search for either "[CX3](=O)[SX1]" or its reverse "[SX1]C(=O)".
      4. On breaking the thioester bond, the fatty acyl fragment must have at least 12 carbons.
      5. The molecule must have at least 4 deprotonated (formally -1 charged) oxygen atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule fits the long-chain fatty acyl-CoA(4-) criteria, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a CoA moiety via an adenine-like fragment.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.GetSubstructMatches(coa_pattern):
        return False, "CoA moiety (adenine substructure) not found"
    
    # 2. Check for thioester bond. 
    # First try pattern where the carbonyl precedes sulfur:
    thioester_pattern1 = Chem.MolFromSmarts("[CX3](=O)[SX1]")
    matches = mol.GetSubstructMatches(thioester_pattern1)
    pattern_used = None
    if matches:
        # Pattern returns a 3-tuple: (carbonyl C, carbonyl O, sulfur)
        pattern_used = "pattern1"
        first_match = matches[0]
        # For pattern1, the acyl chain is attached to the carbonyl carbon.
        acyl_carbon_idx = first_match[0]
        sulfur_idx = first_match[2]
    else:
        # Try reverse pattern: sulfur comes first then carbonyl.
        thioester_pattern2 = Chem.MolFromSmarts("[SX1]C(=O)")
        matches = mol.GetSubstructMatches(thioester_pattern2)
        if matches:
            pattern_used = "pattern2"
            first_match = matches[0]
            # For pattern2, match returns a 2-tuple: (sulfur, carbonyl carbon).
            acyl_carbon_idx = first_match[1]
            sulfur_idx = first_match[0]
        else:
            return False, "No thioester bond (neither [CX3](=O)[SX1] nor [SX1]C(=O)) found"
    
    # 3. Identify the thioester bond by finding the bond between the acyl carbon and the sulfur atom.
    bond = mol.GetBondBetweenAtoms(acyl_carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not correctly identified"
    bond_idx = bond.GetIdx()
    
    # 4. Fragment the molecule by "cutting" the thioester bond.
    frags_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    frags = Chem.GetMolFrags(frags_mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return False, "Could not split molecule into fragments"
    
    # 5. Decide which fragment is the fatty acyl (long-chain) portion.
    # We assume that the CoA moiety is the one that contains the adenine-like substructure.
    fatty_fragment = None
    for frag in frags:
        if not frag.GetSubstructMatches(coa_pattern):
            fatty_fragment = frag
            break
    if fatty_fragment is None:
        # if we cannot discriminate (e.g. both fragments contain adenine or none do), choose the larger fragment
        fatty_fragment = max(frags, key=lambda x: x.GetNumAtoms())
    
    # Count carbons (atomic number 6) in the fatty acyl fragment.
    carbon_count = sum(1 for atom in fatty_fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Fatty acyl chain too short ({carbon_count} carbons found; need at least 12)"
    
    # 6. Check for the deprotonated phosphate groups: at least 4 oxygen atoms with a -1 formal charge overall.
    neg_o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if neg_o_count < 4:
        return False, f"Not enough deprotonated oxygens (found {neg_o_count}; need at least 4 for CoA(4-))"
    
    return True, "Matches long-chain fatty acyl-CoA(4-) criteria"

# Example usage:
if __name__ == "__main__":
    # Test example: palmitoyl-CoA(4-) SMILES
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCCCCCCCCCCCCCC)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
    result, reason = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
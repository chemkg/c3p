"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: Long-chain fatty acyl-CoA(4-)
Definition:
  A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups 
  of any long-chain fatty acyl-CoA; major species at pH 7.3.

Criteria:
  1. The SMILES must be valid.
  2. The molecule must have a CoA moiety as indicated by an adenine-like substructure
     (SMARTS: n1cnc2c(N)ncnc12).
  3. The molecule must contain a thioester bond connecting the fatty acyl chain to the CoA.
     The thioester bond is recognized by searching for either "[CX3](=O)[SX2]" or its
     reverse "[SX2]C(=O)". (Note: sulfur in a thioester is divalent and thus, [SX2] is used.)
  4. Removing the thioester bond should yield a fatty acyl fragment with at least 12 carbons.
  5. The molecule must have at least 4 deprotonated oxygen atoms (O atoms with a formal charge of -1).
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a given molecule (SMILES string) is a long-chain fatty acyl-CoA(4-).

    Criteria:
      1. Valid SMILES.
      2. Contains a CoA moiety detected via an adenine-like fragment (SMARTS: n1cnc2c(N)ncnc12).
      3. Contains a thioester bond connecting the acyl chain and CoA, identified using either
         "[CX3](=O)[SX2]" or its reverse "[SX2]C(=O)" SMARTS.
      4. On severing the thioester bond, the fatty acyl fragment must have at least 12 carbon atoms.
      5. The entire molecule must include at least 4 deprotonated oxygens (O with -1 charge).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule fits the long-chain fatty acyl-CoA(4-) criteria, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a CoA moiety via an adenine-like fragment.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.GetSubstructMatches(coa_pattern):
        return False, "CoA moiety (adenine substructure) not found"
    
    # 2. Look for the thioester bond.
    # Update the SMARTS patterns to match a divalent sulfur in a thioester.
    thioester_pattern1 = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    matches = mol.GetSubstructMatches(thioester_pattern1)
    pattern_used = None
    if matches:
        # For pattern1, we assume a thioester of the type C(=O)-S.
        pattern_used = "pattern1"
        first_match = matches[0]
        # first_match returns (carbonyl carbon, carbonyl oxygen, sulfur)
        acyl_carbon_idx = first_match[0]
        sulfur_idx = first_match[2]
    else:
        # Try the reverse matching: S-C(=O)
        thioester_pattern2 = Chem.MolFromSmarts("[SX2]C(=O)")
        matches = mol.GetSubstructMatches(thioester_pattern2)
        if matches:
            pattern_used = "pattern2"
            first_match = matches[0]
            # For pattern2, the tuple is (sulfur, carbonyl carbon)
            sulfur_idx = first_match[0]
            acyl_carbon_idx = first_match[1]
        else:
            return False, "No thioester bond (neither [CX3](=O)[SX2] nor [SX2]C(=O)) found"
    
    # 3. Identify the thioester bond by finding the bond between acyl carbon and sulfur atom.
    bond = mol.GetBondBetweenAtoms(acyl_carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not correctly identified"
    bond_idx = bond.GetIdx()
    
    # 4. "Cut" the molecule at the thioester bond to separate the fatty acyl fragment.
    frags_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    frags = Chem.GetMolFrags(frags_mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return False, "Could not split molecule into fragments"
    
    # 5. Decide which fragment is the fatty acyl chain.
    # We assume the CoA moiety fragment will contain the adenine-like substructure.
    fatty_fragment = None
    for frag in frags:
        if not frag.GetSubstructMatches(coa_pattern):
            fatty_fragment = frag
            break
    if fatty_fragment is None:
        # If ambiguity exists, choose the largest fragment by atom count.
        fatty_fragment = max(frags, key=lambda x: x.GetNumAtoms())
    
    # Count the number of carbon atoms in the fatty acyl fragment.
    carbon_count = sum(1 for atom in fatty_fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Fatty acyl chain too short ({carbon_count} carbons found; need at least 12)"
    
    # 6. Check that the molecule has at least 4 deprotonated oxygens (oxygen atoms with charge -1).
    neg_o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if neg_o_count < 4:
        return False, f"Not enough deprotonated oxygens (found {neg_o_count}; need at least 4 for CoA(4-))"
    
    return True, "Matches long-chain fatty acyl-CoA(4-) criteria"

# Example usage:
if __name__ == "__main__":
    # Test example: palmitoyl-CoA(4-) SMILES string
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCCCCCCCCCCCCCC)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
    result, reason = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
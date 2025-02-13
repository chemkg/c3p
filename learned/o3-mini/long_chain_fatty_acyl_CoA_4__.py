"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: Long-chain fatty acyl-CoA(4-)
Definition:
  A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups 
  of any long-chain fatty acyl-CoA; major species at pH 7.3.

This classifier checks for:
  - A CoA moiety (using an adenine/inosine-like SMARTS pattern)
  - A thioester bond connecting the acyl chain with CoA. (SMARTS: [CX3](=O)[SX1])
  - A “long-chain” fatty acyl portion; after breaking the thioester bond we count carbons,
    requiring at least 12.
  - At least 4 deprotonated (negatively charged) oxygens in the overall molecule.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a given molecule (SMILES string) is a long-chain fatty acyl-CoA(4-).

    The criteria are:
      1. The molecule must be parsable.
      2. It must contain a CoA moiety; here we look for an adenine-like substructure.
      3. It must contain a thioester bond connecting a fatty acyl chain and CoA,
         using the pattern [CX3](=O)[SX1]. Note that this pattern returns three atoms:
         the carbonyl carbon, the carbonyl oxygen, and the sulfur.
      4. The fatty acyl portion (obtained by cleaving the thioester bond) must have at least 12 carbons.
      5. The molecule must contain at least 4 negatively charged oxygens (deprotonated phosphate groups).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if it fits as long-chain fatty acyl-CoA(4-), False otherwise.
        str: A reason for the classification.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a CoA moiety via an adenine-like fragment.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.GetSubstructMatches(coa_pattern):
        return False, "CoA moiety (adenine substructure) not found"

    # 2. Check for a thioester bond – use [CX3](=O)[SX1].
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX1]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond ([CX3](=O)[SX1]) found"

    # The match returns three atoms: 
    # index 0: thioester carbon (attached to acyl chain),
    # index 1: the carbonyl oxygen,
    # index 2: the sulfur atom.
    first_match = thioester_matches[0]
    acyl_carbon_idx = first_match[0]
    sulfur_idx = first_match[2]

    # 3. Identify the fatty acyl (long-chain) fragment by breaking the thioester bond.
    bond = mol.GetBondBetweenAtoms(acyl_carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not correctly identified"
    
    bond_idx = bond.GetIdx()
    # Fragment the molecule by cutting the thioester bond.
    frags_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    frags = Chem.GetMolFrags(frags_mol, asMols=True, sanitizeFrags=True)
    
    # Decide which fragment is the fatty acyl chain.
    # We choose the fragment that either has a "C(=O)" moiety or is the largest fragment.
    acyl_fragment = None
    for frag in frags:
        frag_smiles = Chem.MolToSmiles(frag)
        if "C(=O)" in frag_smiles or frag_smiles.startswith("C(=O)"):
            acyl_fragment = frag
            break
    if acyl_fragment is None:
        acyl_fragment = max(frags, key=lambda x: x.GetNumAtoms())

    # Count carbon atoms in the acyl fragment.
    carbon_count = sum(1 for atom in acyl_fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Fatty acyl chain too short ({carbon_count} carbons found; need at least 12)"

    # 4. Check for the expected deprotonation of phosphate/diphosphate groups.
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
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
  - A CoA substructure (via the adenine/inosine SMARTS pattern)
  - A thioester bond ([CX3](=O)[SX1]) that connects the acyl chain with the CoA,
  - A “long-chain” acyl portion (we extract the fragment attached to the carbonyl and count its C atoms,
    requiring a minimum count, here set to 12),
  - At least 4 negatively charged oxygen atoms representing deprotonated phosphate groups.
  
If any of these conditions is not met, the function returns False and a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a given molecule (SMILES string) is a long-chain fatty acyl-CoA(4-).

    The criteria are:
      1. The molecule must be parsable.
      2. It must contain a CoA moiety; here we look for the adenine-like substructure.
      3. It must contain a thioester bond connecting a fatty acyl chain and CoA,
         using the pattern [CX3](=O)[SX1].
      4. The fatty acyl portion (the fragment attached at the carbonyl, separate from the CoA)
         must have a “long-chain” (we require at least 12 carbon atoms).
      5. The molecule should show evidence of deprotonation (at least 4 negatively charged oxygens).

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

    # 1. Check for a CoA moiety.
    # Use a simple adenine fragment SMARTS as a proxy for the CoA nucleotide part.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.GetSubstructMatches(coa_pattern):
        return False, "CoA moiety (adenine substructure) not found"

    # 2. Check for thioester bond – a key feature of fatty acyl-CoA.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX1]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond ([CX3](=O)[SX1]) found"

    # For simplicity, we assume one thioester bond. Here we take the first occurrence.
    # The SMARTS match returns a tuple of atom indices: (carbonyl carbon, sulfur)
    acyl_carbon_idx, sulfur_idx = thioester_matches[0]

    # 3. Identify the acyl (fatty) chain fragment.
    # We "break" the bond between the carbonyl carbon and the sulfur.
    bond = mol.GetBondBetweenAtoms(acyl_carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not correctly identified"
    bond_idx = bond.GetIdx()

    # Create a molecule with the thioester bond cut:
    frags_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    # Get individual fragments as separate molecules.
    frags = Chem.GetMolFrags(frags_mol, asMols=True, sanitizeFrags=True)
    
    # We expect that one of these fragments contains the fatty acyl chain.
    # To decide which fragment is the acyl chain, check which contains the acyl carbon.
    acyl_fragment = None
    for frag in frags:
        # If the fragment contains our acyl carbon index, then this is our candidate.
        # Note: FragmentOnBonds creates new atom indices so we match via the original atom mapping.
        # Here we will check if any atom in frag has the same atomic symbol and is in an environment
        # expected for a fatty acyl chain (aliphatic carbon).
        # A simple workaround: generate the SMILES for the fragment and check for aliphatic chain appearance.
        frag_smiles = Chem.MolToSmiles(frag)
        # Our acyl chain should be mostly a string of C’s. Look for a carbonyl at the beginning.
        if frag_smiles.find("C(=O)") != -1 or frag_smiles.startswith("C(=O)"):
            acyl_fragment = frag
            break
    if acyl_fragment is None:
        # As an alternative, take the larger fragment (by number of atoms) which is likely the acyl part.
        acyl_fragment = max(frags, key=lambda x: x.GetNumAtoms())

    # Count carbon atoms in the acyl fragment.
    carbon_count = sum(1 for atom in acyl_fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    # Here we choose a threshold of 12 carbons to call it a "long-chain" fatty acyl
    if carbon_count < 12:
        return False, f"Fatty acyl chain too short ({carbon_count} carbons found; need at least 12)"

    # 4. Check for the expected deprotonation of phosphate/diphosphate groups.
    # Count the number of oxygen atoms with a formal charge of -1.
    neg_o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if neg_o_count < 4:
        return False, f"Not enough deprotonated oxygens (found {neg_o_count}; need at least 4 for CoA(4-))"

    return True, "Matches long-chain fatty acyl-CoA(4-) criteria"
    
# Example usage:
if __name__ == "__main__":
    # Example: palmitoyl-CoA(4-) SMILES
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCCCCCCCCCCCCCC)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
    result, reason = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: Short-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation
of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.
Improved version: Separates the acyl chain by breaking the thioester bond so that
only the acyl (fatty acid) fragment is counted. Also uses several SMARTS patterns 
to identify the CoA moiety.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    The algorithm checks:
      1. That there is at least one thioester bond (C(=O)S).
      2. By breaking the bond between the acyl carbon and sulfur, it isolates the acyl (fatty acid)
         fragment. It then counts the number of carbon atoms (excluding heteroatoms) in this fragment.
         Only fragments with 2 to 6 carbon atoms qualify as “short‐chain”.
      3. Looks for the characteristic adenine substructure (a major part of Coenzyme A) in the full molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Locate the thioester bond: pattern for a carbonyl (C=O) directly attached to an S.
    thioester_smarts = "[C](=O)[S]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) bond detected; not an acyl-CoA"
    
    # We use the first match. The SMARTS is defined so that:
    # match[0] = acyl carbon (the carbon of the carbonyl group),
    # match[1] = the carbonyl oxygen,
    # match[2] = the sulfur.
    acyl_carbon_idx, _, sulfur_idx = thioester_matches[0]
    
    # 2. In order to count the carbons in the acyl chain (the fatty acid fragment),
    # we “cut” the molecule at the thioester bond. This is done by removing the bond from
    # the acyl carbon to the sulfur so that the fatty acyl fragment separates from the CoA.
    rw_mol = Chem.RWMol(mol)
    bond = rw_mol.GetBondBetweenAtoms(acyl_carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Unexpected error: thioester bond not found in molecule"
    rw_mol.RemoveBond(acyl_carbon_idx, sulfur_idx)
    # Get the fragments (each is given as a tuple of atom indices from the original molecule).
    frags = Chem.GetMolFrags(rw_mol, asMols=False)
    
    # Identify the fragment that contains the acyl carbon.
    acyl_frag = None
    for frag in frags:
        if acyl_carbon_idx in frag:
            acyl_frag = frag
            break
    if acyl_frag is None:
        return False, "Error in fragmenting molecule by breaking thioester bond"
    
    # Count the number of carbon atoms (atomic number 6) in the acyl fragment.
    acyl_carbons = sum(1 for idx in acyl_frag if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if acyl_carbons < 2:
        return False, f"Acyl fragment has only {acyl_carbons} carbon(s); must be at least 2 carbons long"
    if acyl_carbons > 6:
        return False, f"Acyl fragment has {acyl_carbons} carbons; too long for a short-chain fatty acid (max 6 carbons)"
    
    # 3. Check for key structural elements of Coenzyme A.
    # One strong indicator is the adenine ring. We test several SMARTS variations.
    adenine_smarts_list = [
        "c1ncnc2ncnc12",       # generic purine pattern (adenine)
        "n1cnc2nc[nH]c2n1",
        "c1nc2nc[nH]c2n1"
    ]
    adenine_found = False
    for smarts in adenine_smarts_list:
        adenine_pat = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(adenine_pat, useChirality=False):
            adenine_found = True
            break
    if not adenine_found:
        return False, "No adenine (CoA) moiety detected"
    
    return True, f"Found thioester bond with an acyl fragment of {acyl_carbons} carbons and a CoA (adenine) moiety"

# Optional: Testing examples
if __name__ == "__main__":
    # Test one example: (R)-3-hydroxypentanoyl-CoA (acyl chain with 5 carbons)
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C[C@@H](CC)O)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
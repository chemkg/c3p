"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA

A medium-chain fatty acyl-CoA is defined as a fatty acyl-CoA that results from the
formal condensation of the thiol group of coenzyme A with the carboxy group of any
medium-chain fatty acid.
Heuristic requirements:
  (i) the molecule contains a CoA-like fragment,
  (ii) a thioester linkage [C(=O)S] linking the fatty acyl chain to CoA,
  (iii) after “cutting” the thioester bond, the isolated acyl fragment contains
       between 6 and 12 carbon atoms and is aliphatic (i.e. no rings or aromatic carbons).
       
Note: This is a rough heuristic that may reject some edge cases.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium‐chain fatty acyl‐CoA based on its SMILES.
    
    Checks for:
      - A CoA-like fragment (using a SMARTS that appears in many acyl-CoA molecules).
      - A thioester group ([CX3](=O)[S]) linking the acyl chain to CoA.
      - The acyl fragment, obtained by “fragmenting” the molecule at the thioester bond,
        should contain between 6 and 12 carbon atoms and no rings or aromatic carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if classified as medium-chain fatty acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a CoA-like fragment.
    # This SMARTS is a rough fragment seen in many acyl-CoA molecules.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A fragment not detected"
    
    # Look for the thioester group linking the acyl chain and CoA.
    thioester_smarts = "[CX3](=O)[S]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    th_matches = mol.GetSubstructMatches(thioester_pattern)
    if not th_matches:
        return False, "Thioester group (C(=O)S) not found"
    
    # Take the first thioester match.
    # The match gives the indices of the carbonyl carbon, the oxygen, and the sulfur.
    carbonyl_idx, _, sulfur_idx = th_matches[0]
    
    # Mark the carbonyl atom so we can track which fragment contains the acyl chain.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    carbonyl_atom.SetProp("acyl_marker", "yes")
    
    # Identify the bond between the carbonyl carbon and the sulfur.
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not properly defined"
    bond_idx = bond.GetIdx()
    
    # Fragment the molecule on this bond, cutting the acyl chain off from the CoA.
    fragmented = Chem.FragmentOnBonds(mol, [bond_idx])
    frags = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
    
    # Identify the acyl fragment as the one that contains our marked carbonyl atom.
    acyl_frag = None
    for frag in frags:
        for atom in frag.GetAtoms():
            if atom.HasProp("acyl_marker"):
                acyl_frag = frag
                break
        if acyl_frag is not None:
            break
    if acyl_frag is None:
        return False, "Could not isolate acyl (fatty acid) fragment from the thioester bond"
    
    # Count carbon atoms in the acyl fragment (ignore dummy atoms with atomic number 0).
    carbon_count = sum(1 for atom in acyl_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Heuristic: in fatty acid chains the carbonyl carbon is part of the chain.
    if not (6 <= carbon_count <= 12):
        return False, f"Acyl chain has {carbon_count} carbons, not in the medium-chain range (6-12)"
    
    # Additional check: fatty acyl chains are typically linear and aliphatic.
    # Reject if any carbon atom in the acyl fragment is aromatic or part of a ring.
    for atom in acyl_frag.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if atom.GetIsAromatic():
                return False, "Acyl fragment contains aromatic carbon(s), not a typical fatty acyl chain"
            if atom.IsInRing():
                return False, "Acyl fragment contains ring atoms, not a typical linear fatty acyl chain"
    
    return True, f"Detected medium-chain fatty acyl-CoA (acyl chain with {carbon_count} carbons)"

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one example: hexanoyl-CoA
    test_smiles = "CCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, explanation = is_medium_chain_fatty_acyl_CoA(test_smiles)
    print(result, explanation)
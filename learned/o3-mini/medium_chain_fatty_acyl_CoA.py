"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA

A medium-chain fatty acyl-CoA is defined as a fatty acyl-CoA that results from the 
formal condensation of the thiol group of coenzyme A with the carboxy group of any medium‐chain fatty acid.
We heuristically require that:
  (i) the molecule contains a CoA-like fragment,
  (ii) a thioester linkage [C(=O)S] linking the fatty acyl chain to CoA,
  (iii) after “cutting” the thioester bond, the acyl fragment (the fatty acid part) contains between 6 and 12 carbon atoms.
  
Note: This annotation is a rough heuristic and may not capture every edge case.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium‐chain fatty acyl‐CoA based on its SMILES string.
    
    The algorithm checks for:
      - A CoA fragment (using a SMARTS pattern that matches a characteristic substructure).
      - A thioester group (pattern "[CX3](=O)[S]") linking the fatty acyl chain to CoA.
      - The acyl fragment isolated from fragmenting at the thioester bond is then analyzed for the number of carbon atoms.
        We consider a fatty acyl chain medium‐chain if it contains between 6 and 12 carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as medium‐chain fatty acyl‐CoA, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First check for a CoA fragment.
    # Here we use a SMARTS that is found in many acyl-CoA molecules: the pattern "SCCNC(=O)CCNC(=O)"
    # Note: This pattern does not capture all possible CoA variations but provides a useful heuristic.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A fragment not detected"
    
    # Look for the thioester group linking the acyl chain and the CoA.
    # This pattern matches a carbonyl carbon attached to a sulfur: C(=O)S
    thioester_smarts = "[CX3](=O)[S]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    th_matches = mol.GetSubstructMatches(thioester_pattern)
    if not th_matches:
        return False, "Thioester group (C(=O)S) not found"

    # Use the first thioester match.
    # th_matches returns tuples of atom indices corresponding to [C, O, S] in the pattern.
    carbonyl_idx, _, sulfur_idx = th_matches[0]

    # Mark the carbonyl atom (i.e. the acyl carbon) so that we can track it after fragmentation.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    carbonyl_atom.SetProp("acyl_marker", "yes")
    
    # Identify the bond between the carbonyl carbon and the sulfur.
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not properly defined"
    bond_idx = bond.GetIdx()
    
    # Fragment the molecule on the thioester bond.
    # This will “cut” the molecule between the fatty acyl chain and the CoA moiety.
    fragmented = Chem.FragmentOnBonds(mol, [bond_idx])
    frags = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
    
    # Identify the acyl fragment by finding the fragment that contains our marked carbonyl atom.
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
    
    # Count carbon atoms in the acyl fragment.
    # (Ignore dummy atoms that result from bond fragmentation; these have atomic number 0.)
    carbon_count = 0
    for atom in acyl_frag.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
            
    # In many fatty acids the carbonyl carbon is counted as part of the chain.
    # We define medium-chain roughly as having between 6 and 12 carbons.
    if not (6 <= carbon_count <= 12):
        return False, f"Acyl chain has {carbon_count} carbons, not in the medium-chain range (6-12)"
    
    return True, f"Detected medium-chain fatty acyl-CoA (acyl chain with {carbon_count} carbons)"

# Example usage:
if __name__ == "__main__":
    test_smiles = "CCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"  # hexanoyl-CoA
    is_acyl_coa, reason = is_medium_chain_fatty_acyl_CoA(test_smiles)
    print(is_acyl_coa, reason)
"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA

A medium-chain fatty acyl-CoA is defined as a fatty acyl-CoA that results from the
formal condensation of the thiol group of coenzyme A with the carboxy group of any
medium-chain fatty acid.

Heuristic checks:
  (i) The molecule must contain a CoA-like fragment; here we require detection of
      the adenine substructure (a decorated purine ring) commonly present in CoA.
  (ii) A thioester linkage [C(=O)S] must be present.
  (iii) After “cutting” the thioester bond, the isolated acyl fragment (the one that
        contains the carbonyl but not the sulfur from the thioester) should be aliphatic,
        with no aromatic atoms or rings, and should have a total of 6–12 carbon atoms.
        
Note: This heuristic is rough and may reject borderline cases.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium‐chain fatty acyl‐CoA based on its SMILES.

    Steps:
      1. Parse the molecule.
      2. Check for a CoA-like fragment by requiring the adenine substructure.
         A loose SMARTS "c1nc2c(n1)nc(nc2)" is used to catch decorated adenine rings.
      3. Identify the thioester bond using the SMARTS pattern [CX3](=O)[S].
      4. Mark the carbonyl atom (acyl side) and the sulfur atom (CoA side).
      5. Fragment the molecule on the thioester bond.
      6. Identify the fragment containing the acyl (fatty acid) chain: it should contain
         the carbonyl marker but not the thioester sulfur marker.
      7. Count the number of carbon atoms in that fragment.
      8. Verify that the acyl fragment is aliphatic (i.e. without aromatic carbons or rings).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if classified as a medium-chain fatty acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    
    # Step 1: Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2: Look for a CoA-like fragment by detecting an adenine substructure.
    # Use a less strict SMARTS so variations like "n1cnc2c(N)ncnc12" are detected.
    adenine_smarts = "c1nc2c(n1)nc(nc2)"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if adenine is None or not mol.HasSubstructMatch(adenine):
        return False, "Coenzyme A adenine substructure not detected"
    
    # Step 3: Locate the thioester group [CX3](=O)[S]
    thioester_smarts = "[CX3](=O)[S]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    th_matches = mol.GetSubstructMatches(thioester)
    if not th_matches:
        return False, "Thioester group (C(=O)S) not found"
    
    # For this heuristic, we take the first thioester match.
    # Each match yields a tuple of three atom indices: (carbonyl carbon, oxygen of C=O, sulfur)
    carbonyl_idx, _, sulfur_idx = th_matches[0]
    
    # Step 4: Mark the key atoms in the molecule for later identification.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
    carbonyl_atom.SetProp("acyl_marker", "yes")
    sulfur_atom.SetProp("thioester_S", "yes")
    
    # Identify the bond connecting the carbonyl carbon and the sulfur atom.
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not properly defined"
    bond_idx = bond.GetIdx()
    
    # Step 5: Fragment the molecule on the thioester bond.
    fragmented = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    frags = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
    
    # Step 6: Determine which fragment corresponds to the acyl (fatty acid) chain.
    # We choose the fragment that has the carbonyl marker (acyl side) but not the sulfur marker.
    acyl_frag = None
    for frag in frags:
        has_acyl_marker = any(atom.HasProp("acyl_marker") for atom in frag.GetAtoms())
        has_thioester_S = any(atom.HasProp("thioester_S") for atom in frag.GetAtoms())
        if has_acyl_marker and not has_thioester_S:
            acyl_frag = frag
            break
    if acyl_frag is None:
        return False, "Could not isolate acyl (fatty acid) fragment from the thioester bond"
    
    # Step 7: Count carbon atoms in the acyl fragment.
    carbon_count = sum(1 for atom in acyl_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (6 <= carbon_count <= 12):
        return False, f"Acyl chain has {carbon_count} carbons, not in the medium-chain range (6-12)"
    
    # Step 8: Ensure that the acyl fragment is aliphatic:
    # It should not contain aromatic carbons or any atoms in rings.
    for atom in acyl_frag.GetAtoms():
        if atom.GetAtomicNum() == 6:  # check carbons
            if atom.GetIsAromatic():
                return False, "Acyl fragment contains aromatic carbon(s), not typical for fatty acyl chains"
            if atom.IsInRing():
                return False, "Acyl fragment contains ring atoms, expected a linear (or branched) chain"
    
    return True, f"Detected medium-chain fatty acyl-CoA with an acyl chain of {carbon_count} carbons"

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one example: hexanoyl-CoA.
    test_smiles = "CCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, explanation = is_medium_chain_fatty_acyl_CoA(test_smiles)
    print(result, explanation)
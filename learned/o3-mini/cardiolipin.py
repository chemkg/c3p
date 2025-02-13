"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: Cardiolipin 
Definition: A phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a molecule of glycerol.
This heuristic method checks:
  • Exactly 2 phosphorus atoms,
  • At least 4 fatty acid (acyl) ester groups (using a SMARTS pattern),
  • The presence of a glycerol linker (HOCH2–CHOH–CH2OH motif) that is expected in cardiolipin,
  • A high molecular weight (typically >1000 Da),
  • And a sufficiently high number of rotatable bonds.
Note: This classification is heuristic and may not cover all examples.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    Cardiolipin is defined as a phosphatidylglycerol composed of two phosphatidic acids 
    linked via a glycerol molecule (hence, usually containing 2 P atoms, 4 acyl chains, 
    and a characteristic glycerol bridge).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as cardiolipin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule has exactly 2 phosphorus (P, atomic number = 15) atoms.
    phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphorus_count != 2:
        return False, f"Expected 2 phosphorus atoms, but found {phosphorus_count}"
    
    # Check for at least 4 acyl ester groups.
    # The SMARTS "OC(=O)CCCC" looks for an ester function with a fatty acid chain (at least 4 carbons).
    acyl_chain_pattern = Chem.MolFromSmarts("OC(=O)CCCC")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) < 4:
        return False, f"Expected at least 4 acyl chain ester groups, but found {len(acyl_matches)}"
    
    # Check for high molecular weight. Cardiolipins are large lipids (typically >1000 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low for cardiolipin: {mol_wt:.2f} Da"
    
    # Check that the molecule has a sufficient number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Not enough rotatable bonds ({n_rotatable}) expected for a flexible cardiolipin"
    
    # Additional critical check: verify the presence of a central glycerol linker.
    # Cardiolipin is defined by a glycerol unit bridging two phosphatidic acid moieties.
    # Here we search for a fragment resembling glycerol: HOCH2-CHOH-CH2OH.
    # The SMARTS below “OCC(O)CO” is a simple approximation for such a motif.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol bridging moiety not found, which is required for cardiolipin"
    
    # If all tests pass, we classify the molecule as cardiolipin.
    return True, "Molecule matches the heuristic criteria for cardiolipin (2 P atoms, ≥4 acyl chains, glycerol linker, high molecular weight)"
    
# Example usage (if this module is run directly):
if __name__ == "__main__":
    # Test with one of the provided cardiolipin examples
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(OC[C@@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)=O)(O)=O"
    result, reason = is_cardiolipin(test_smiles)
    print(f"Result: {result}\nReason: {reason}")
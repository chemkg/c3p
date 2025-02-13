"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: Cardiolipin (A phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a molecule of glycerol)
This heuristic method checks for:
  • exactly 2 phosphorus atoms,
  • at least 4 acyl chain ester groups (heuristically defined using a SMARTS pattern),
  • a sufficiently high molecular weight (typical cardiolipins are large molecules, e.g. >1000 Da).
Note: This classification is heuristic and may not cover all examples.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    Cardiolipin is defined as a phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a glycerol.
    In practice, cardiolipin molecules contain two phosphate groups and four fatty acid (acyl) chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule passes heuristic cardiolipin criteria, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule contains exactly 2 phosphorus atoms (P atomic number = 15)
    phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphorus_count != 2:
        return False, f"Expected 2 phosphorus atoms, but found {phosphorus_count}"
    
    # Heuristic check for fatty acid (acyl chain) ester groups.
    # We use a SMARTS that looks for an ester motif with a carbonyl connected to a chain of at least 4 carbons.
    # Note: This is a simplification; cardiolipin molecules contain 4 acyl chains.
    acyl_chain_pattern = Chem.MolFromSmarts("OC(=O)CCCC")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) < 4:
        return False, f"Expected at least 4 acyl chain ester groups, but found {len(acyl_matches)}"
    
    # Check that the molecular weight is high enough (cardiolipins are large lipids, typically above 1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low for cardiolipin: {mol_wt:.2f} Da"
    
    # Additional checks (e.g., rotatable bonds) could be added here if desired
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Not enough rotatable bonds ({n_rotatable}) expected for a flexible cardiolipin"
    
    return True, "Molecule matches the heuristic criteria for cardiolipin (2 P atoms, ≥4 acyl chains, high molecular weight)"

# Example usage (if running this script directly):
if __name__ == "__main__":
    # Example SMILES for a cardiolipin (one of the provided examples)
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(OC[C@@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)=O)(O)=O"
    result, reason = is_cardiolipin(test_smiles)
    print(f"Result: {result}\nReason: {reason}")
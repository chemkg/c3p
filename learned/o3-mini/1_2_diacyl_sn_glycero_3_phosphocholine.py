"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
#!/usr/bin/env python3
"""
Classifies: 1,2-diacyl-sn-glycero-3-phosphocholine
Definition: The conjugate base of a 1,2-diacyl-sn-glycero-3-phosphocholine compound formed by deprotonation of the phosphate OH group.
A typical structure of these molecules (phosphatidylcholines) is a glycerol backbone with fatty acyl (ester) chains at the sn-1 and sn-2 positions and a phosphocholine headgroup at the sn-3 position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a 1,2-diacyl-sn-glycero-3-phosphocholine.
    
    This function checks two main requirements:
      1. It finds at least two ester bonds (OC(=O)) that are indicative of two fatty acyl chains.
      2. It finds a phosphocholine headgroup, that is, a phosphate group bearing a negative charge and attached to a choline fragment (OCC[N+](C)(C)C).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is recognized as a 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a phosphocholine headgroup.
    # We look for a phosphate (P=O with a deprotonated oxygen, [O-]) attached to an oxygen and then 
    # to a choline fragment "OCC[N+](C)(C)C". 
    phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine headgroup not found"
    
    # Check for the presence of at least two ester groups (indicative of two acyl chains).
    # Ester bonds are typically represented as OC(=O) - an oxygen attached to a carbonyl.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} acyl ester(s); at least 2 are required"
    
    # Optional sanity check: verify that there is at least one phosphorus atom present.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) == 0:
        return False, "No phosphorus atom found"
    
    # Additional optional checks include verifying the presence of a choline substructure.
    choline_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if not mol.HasSubstructMatch(choline_pattern):
        return False, "Choline moiety not found"
    
    # If all required patterns are found, return True.
    return True, "Molecule contains a phosphocholine headgroup with two acyl ester chains consistent with a 1,2-diacyl-sn-glycero-3-phosphocholine"
    
# Example usage (uncomment for testing):
# test_smiles = "P(OC[C@@H](COC(CCCCCCCCCCCCCCCCC)=O)OC(=O)CCCC)(=O)(OCC[N+](C)(C)C)[O-]"
# is_pc, reason = is_1_2_diacyl_sn_glycero_3_phosphocholine(test_smiles)
# print(is_pc, reason)
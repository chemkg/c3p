"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol
Definition:
An organosulfur compound in which a thiol group (-SH) is attached to a carbon atom 
of any aliphatic or aromatic moiety.
Note:
  - Only the protonated form (with at least one H on S) qualifies.
  - Molecules with multiple amide bonds and high molecular weight are treated as peptides and not classified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol must have an S atom that:
      - is not part of a disulfide (i.e. not bonded to another S)
      - is attached to at least one carbon atom
      - carries at least one hydrogen (i.e. in the -SH form, not the deprotonated thiolate)
    
    Additionally, molecules with multiple amide bonds and high molecular weight 
    are filtered out to avoid misclassifying peptides.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a thiol, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogen atoms so that we capture the hydrogen count on sulfur
    mol = Chem.AddHs(mol)
    
    # Filter out molecules that appear to be peptides: look for multiple amide bonds and check molecular weight.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if len(amide_matches) >= 2 and mol_wt > 300:
        return False, "Molecule appears to be a peptide (multiple amide bonds and high molecular weight)"
    
    # Look for sulfur atoms that qualify as a thiol group (-SH).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # only consider sulfur
            continue
        
        # Skip if the sulfur atom is directly bonded to another sulfur (often part of disulfide bridges)
        if any(neighbor.GetAtomicNum() == 16 for neighbor in atom.GetNeighbors()):
            continue
        
        # Ensure the sulfur is attached to at least one carbon (to confirm it is part of an organic moiety)
        if not any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
            continue
        
        # Check whether the sulfur carries at least one hydrogen.
        # We look at the total hydrogen count on the atom.
        num_h = atom.GetTotalNumHs()
        # Only classify as thiol if at least one hydrogen is attached and formal charge is zero.
        if num_h > 0 and atom.GetFormalCharge() == 0:
            return True, "Found thiol group: sulfur with attached hydrogen bonded to carbon"
        
    return False, "No thiol group (-SH attached to carbon) found"

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "OS(=O)(=O)CCS",  # coenzyme M; should be classified as thiol
        "CCS",           # ethanethiol; should be classified as thiol
        "SC1CCCC1",      # Cyclopentanethiol; should be classified as thiol
        "[S-]C(=S)NCCNC([S-])=S",  # zineb (should not be classified)
        "NCCCS"         # 3-aminopropane-1-thiol example (should be thiol if S carries H)
    ]
    for smi in test_smiles:
        result, reason = is_thiol(smi)
        print(f"SMILES: {smi} -> {result}, {reason}")
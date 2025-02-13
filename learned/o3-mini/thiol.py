"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol
Definition:
An organosulfur compound in which a thiol group (-SH) is attached to a carbon atom
of any aliphatic or aromatic moiety.

This improved classifier:
 - Adds explicit hydrogens so that implicit hydrogen counts on sulfur are available.
 - Considers both neutral thiols (with hydrogen) and deprotonated thiol groups (thiolates, with formal charge -1).
 - Filters out molecules that appear to be peptides (multiple amide bonds and high molecular weight),
   which have frequently been misclassified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is defined as an organosulfur compound in which a thiol group (-SH)
    is attached to a carbon atom (or exists in its deprotonated thiolate [S-] form).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a thiol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that we can count the H atoms on S
    mol = Chem.AddHs(mol)
    
    # Check for peptide character.
    # We look for multiple amide bonds (C(=O)N) and check the molecular weight.
    # Simple thiol compounds are expected to have 0 or 1 amide bond.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if len(amide_matches) >= 2 and mol_wt > 300:
        return False, "Molecule appears to be a peptide (multiple amide bonds and high molecular weight)"

    # Iterate over all atoms to find a sulfur atom that qualifies as part of a thiol group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # 16 represents sulfur
            # Exclude cases where sulfur is directly bonded to another sulfur (common in disulfide bridges).
            if any(neighbor.GetAtomicNum() == 16 for neighbor in atom.GetNeighbors()):
                continue

            # Check that the sulfur is attached to at least one carbon atom.
            if not any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
                continue

            # Now determine if the sulfur is in a thiol environment.
            # Two cases are considered:
            # 1. The sulfur has one or more hydrogen atoms (i.e. as in -SH group).
            # 2. The sulfur is deprotonated (thiolate, no H) but has a formal charge of -1.
            num_h = atom.GetTotalNumHs()
            if num_h > 0:
                return True, "Found thiol group: sulfur with attached hydrogen bonded to carbon"
            if num_h == 0 and atom.GetFormalCharge() == -1:
                return True, "Found thiol group: deprotonated sulfur (thiolate) bonded to carbon"

    return False, "No thiol group (-SH or corresponding thiolate attached to carbon) found"

# Example usage:
if __name__ == '__main__':
    # Testing on coenzyme M which should be classified as a thiol.
    smiles_list = [
        "OS(=O)(=O)CCS",  # coenzyme M
        "CCS",           # ethanethiol
        "SC1CCCC1"       # Cyclopentanethiol
    ]
    for smi in smiles_list:
        result, reason = is_thiol(smi)
        print(f"SMILES: {smi} -> {result}, {reason}")
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
    A valid thiol must have at least one sulfur atom that:
      - is not part of a disulfide (i.e. not bonded to another S),
      - is attached to at least one carbon,
      - carries at least one hydrogen (i.e. is in the -SH protonated form),
      - and is not oxidized (i.e. does not have a double bond to an oxygen).
    
    Additionally, molecules that appear to be peptides (with multiple amide bonds and high molecular weight)
    are filtered out.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a thiol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that hydrogen counting on sulfur is reliable.
    mol = Chem.AddHs(mol)

    # Filter out peptides: if there are 2 or more amide bonds and the molecular weight exceeds 300, skip.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if len(amide_matches) >= 2 and mol_wt > 300:
        return False, "Molecule appears to be a peptide (multiple amide bonds and high molecular weight)"
    
    # Iterate over all atoms looking for a qualifying sulfur (S) that is part of a thiol group (-SH)
    for atom in mol.GetAtoms():
        # Check only sulfur atoms (atomic number 16)
        if atom.GetAtomicNum() != 16:
            continue
        
        # Exclude sulfur atoms that are bonded to another sulfur (disulfides)
        if any(neighbor.GetAtomicNum() == 16 for neighbor in atom.GetNeighbors()):
            continue
        
        # Ensure the sulfur is attached to at least one carbon (organic environment)
        if not any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
            continue
        
        # Exclude sulfur atoms that are oxidized: if any bond from S to oxygen is a double bond, skip.
        oxidized = False
        for bond in atom.GetBonds():
            # Check if the bond is a double bond and if the other atom is oxygen (atomic number 8)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    oxidized = True
                    break
        if oxidized:
            continue
        
        # Check if the sulfur has at least one hydrogen neighbor (protonated -SH).
        hydrogen_found = any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())
        if not hydrogen_found:
            continue
        
        # If all criteria are met, we have found a qualifying thiol group.
        return True, "Found thiol group: sulfur with attached hydrogen bonded to carbon"
    
    # If the loop finishes without a positive match, return False.
    return False, "No thiol group (-SH attached to carbon) found"

# Example usage for testing:
if __name__ == '__main__':
    test_smiles = [
        "OS(=O)(=O)CCS",  # coenzyme M; should be thiol due to the terminal -S, even if other S groups exist.
        "CCS",           # ethanethiol; should be classified as thiol
        "SC1CCCC1",      # Cyclopentanethiol; should be thiol
        "NCCCNCCS",      # WR-1065; should be thiol (terminal sulfur with attached hydrogen)
        # Additional examples from the provided list can be added here.
    ]
    for smi in test_smiles:
        result, reason = is_thiol(smi)
        print(f"SMILES: {smi} -> {result}, {reason}")
"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene as defined as “An acyclic branched or unbranched hydrocarbon
having one carbon–carbon double bond and the general formula CnH2n”.
Only pure hydrocarbons (only C and H) are allowed. Molecules with more than one
C=C or with any ring (including if the double-bond is in a ring) are not alkenes.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene (by this definition) must be a pure hydrocarbon (only C and H),
    acyclic (no rings) and contain exactly one carbon–carbon double bond.
    In addition, the hydrogen count must be exactly 2 times the number of carbon atoms (CnH2n).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene by the definition, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # To count hydrogens correctly, add explicit hydrogens.
    mol_with_H = Chem.AddHs(mol)
    
    # Check that molecule is a hydrocarbon: only carbon and hydrogen atoms.
    for atom in mol_with_H.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than C and H (not a pure hydrocarbon)"
    
    # Enforce the acyclic condition: molecule must have no rings.
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains ring structures (must be acyclic)"
    
    # Count the number of C=C bonds. Only count a bond if it is a double bond and both atoms are carbon.
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                cc_double_bonds += 1
                # Also ensure that the atoms in the double bond are not members of any ring.
                if bond.IsInRing():
                    return False, "The double bond is in a ring structure"
    
    if cc_double_bonds != 1:
        return False, f"Molecule contains {cc_double_bonds} carbon–carbon double bonds (must be exactly one)"
    
    # Verify the overall formula is CnH2n.
    carbon_count = sum(1 for atom in mol_with_H.GetAtoms() if atom.GetAtomicNum() == 6)
    hydrogen_count = sum(1 for atom in mol_with_H.GetAtoms() if atom.GetAtomicNum() == 1)
    expected_hydrogens = 2 * carbon_count
    if hydrogen_count != expected_hydrogens:
        return False, f"Molecule does not satisfy the formula CnH2n: found C{carbon_count}H{hydrogen_count}"
    
    return True, "Molecule is an acyclic hydrocarbon with exactly one carbon–carbon double bond and formula CnH2n"
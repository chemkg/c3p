"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene as defined as “An acyclic branched or unbranched hydrocarbon
having one carbon–carbon double bond and the general formula CnH2n”.
Only pure hydrocarbons (only C and H) are allowed. Molecules with more than one
C=C (or with any ring, including if the double-bond is in a ring) are not alkenes.
"""

import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene (by this definition) must be a pure hydrocarbon (only C and H),
    acyclic (no rings) and contain exactly one carbon–carbon double bond.
    In addition, the molecular formula must be exactly CnH2n.
    
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
        
    # Check that the molecule is a pure hydrocarbon: only carbon (6) and hydrogen (1) atoms.
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (1, 6):
            return False, "Molecule contains atoms other than C and H (not a pure hydrocarbon)"
    
    # Enforce the acyclic condition: molecule must have no rings.
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains ring structures (must be acyclic)"
        
    # Count the number of carbon-carbon double bonds. 
    # Walk through all bonds; if the bond is a double bond and is between two carbons, count it.
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Both ends must be carbon.
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                # Even if the molecule overall is acyclic, double-check that the double bond is not part of any ring.
                if bond.IsInRing():
                    return False, "The carbon–carbon double bond is in a ring structure (must be acyclic)"
                cc_double_bonds += 1
                
    if cc_double_bonds != 1:
        return False, f"Molecule contains {cc_double_bonds} carbon–carbon double bonds (must be exactly one)"
        
    # Verify that the overall molecular formula is exactly CnH2n.
    # We use rdMolDescriptors.CalcMolFormula which returns a string like 'C10H20'
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # The expected formula should have only carbon and hydrogen.
    match = re.fullmatch(r"C(\d*)(H)(\d+)", formula)
    if not match:
        return False, f"Molecular formula {formula} does not match expected pattern CnH2n"
    
    # Get the carbon count; if no number appears after C, it means 1.
    c_count_str = match.group(1)
    c_count = int(c_count_str) if c_count_str != "" else 1
    h_count = int(match.group(3))
    
    # Check if hydrogen count equals 2 * carbon count.
    if h_count != 2 * c_count:
        return False, f"Molecule does not satisfy the formula CnH2n: found C{c_count}H{h_count}"
    
    return True, "Molecule is an acyclic hydrocarbon with exactly one carbon–carbon double bond and formula CnH2n"
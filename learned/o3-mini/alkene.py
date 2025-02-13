"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene

An alkene (by this definition) is an acyclic (ring‐free) pure hydrocarbon
(i.e. only carbon and hydrogen) that contains exactly one carbon–carbon double bond.
In addition, its overall formula must be CnH2n, which is only possible if the hydrogen
deficiency comes exclusively from that one C=C.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    Criteria:
      • Pure hydrocarbon: contains only carbon (atomic number 6) and hydrogen (atomic number 1).
      • Acyclic: no rings are present.
      • Contains exactly one carbon–carbon double bond. (Double bonds in rings are not allowed.)
      • The overall formula must be CnH2n.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule meets the alkene definition, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemical annotations to avoid ambiguities.
    Chem.RemoveStereochemistry(mol)
    
    # Check that every atom is either carbon (6) or hydrogen (1).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than C and H (not a pure hydrocarbon)"
    
    # Enforce acyclicity: reject if any ring is found.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures (must be acyclic)"
    
    # Count the number of carbon–carbon double bonds.
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        # We only consider double bonds between two carbons.
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Even though the molecule is overall acyclic, double-check that this bond is not in a ring.
                if bond.IsInRing():
                    return False, "A carbon–carbon double bond is in a ring (must be acyclic)"
                cc_double_bonds += 1
    if cc_double_bonds != 1:
        return False, f"Molecule contains {cc_double_bonds} carbon–carbon double bond(s) (must be exactly one)"
    
    # To verify the overall formula is CnH2n, we add explicit hydrogens and count atoms.
    mol_H = Chem.AddHs(mol)
    c_count = sum(1 for atom in mol_H.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol_H.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if h_count != 2 * c_count:
        return False, f"Molecule does not satisfy the formula CnH2n: found C{c_count}H{h_count}"
    
    return True, "Molecule is an acyclic pure hydrocarbon with exactly one carbon–carbon double bond and formula CnH2n"


# Example usage:
if __name__ == "__main__":
    # Test with a few examples
    examples = [
        ("CCCC\\C=C\\C", "(E)-2-octene"),
        ("CCCCCCCCCCCCCCC\\C=C\\CCCCCC", "(7E)-tricosene"),
        ("C=CCCCCCCC1OC1", "1,2-epoxy-9-decene"),
        ("C(CCCC(CCCC(C=C)=C)C)(CCCC(C)C)C", "neophytadiene"),
    ]
    
    for smi, name in examples:
        result, reason = is_alkene(smi)
        print(f"{name} (SMILES: {smi}) -> {result}. Reason: {reason}")
"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene

An alkene (by this definition) is an acyclic (ring‐free) pure hydrocarbon 
(i.e. only carbon and hydrogen) that contains exactly one carbon–carbon double bond. 
In addition, its overall formula must be CnH2n (i.e. the only unsaturation comes 
from that one double bond).
"""

import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    Criteria:
      • Pure hydrocarbon: contains only carbon (atomic number 6) and hydrogen (atomic number 1)
      • Acyclic: no rings are present.
      • Contains exactly one carbon–carbon double bond (and that bond is not in a ring).
      • The overall molecular formula (including implicit hydrogens) must be CnH2n.
        (This condition is equivalent to having exactly one degree of unsaturation 
         for an acyclic hydrocarbon.)
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule meets the alkene definition, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove stereochemical annotations to avoid ambiguities.
    Chem.RemoveStereochemistry(mol)

    # Ensure that every atom is either carbon (atomic number 6) or hydrogen (atomic number 1)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than C and H (not a pure hydrocarbon)"

    # Reject molecules that contain any rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures (must be acyclic)"

    # Count the number of carbon–carbon double bonds.
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Even though we already checked acyclicity, double-check that the double bond is not in a ring.
                if bond.IsInRing():
                    return False, "A carbon–carbon double bond is in a ring (must be acyclic)"
                cc_double_bonds += 1

    if cc_double_bonds != 1:
        return False, f"Molecule contains {cc_double_bonds} carbon–carbon double bond(s) (must be exactly one)"

    # Calculate the molecular formula using RDKit.
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Expect formula of the form C\d+H\d+ (e.g. "C8H16")
    m = re.fullmatch(r"C(\d+)H(\d+)", formula)
    if not m:
        return False, f"Unexpected molecular formula format: {formula}"
    c_count = int(m.group(1))
    h_count = int(m.group(2))

    # The condition for the alkene is that the formula must be CnH2n.
    if h_count != 2 * c_count:
        return False, f"Molecule does not satisfy the formula CnH2n: found {formula}"

    # Alternatively, compute degree of unsaturation.
    # For an acyclic hydrocarbon, unsaturation = (2C + 2 - H) / 2.
    unsat = (2 * c_count + 2 - h_count) / 2
    if unsat != 1:
        return False, f"Molecule has {unsat} degree(s) of unsaturation (expected exactly one for an alkene)"

    return True, "Molecule is an acyclic pure hydrocarbon with exactly one C=C double bond and formula CnH2n"

# Example usage (you can run tests when executing this file directly):
if __name__ == "__main__":
    test_examples = [
        ("CCCC\\C=C\\C", "(E)-2-octene"),
        ("CCCCCCCCCCCCCCC\\C=C\\CCCCCC", "(7E)-tricosene"),
        ("CCC(C)CCCCCCCCCCCCCCC=C", "17-methylnonadec-1-ene"),
        ("[H]C(C)=C([H])CCCCC", "2-octene"),
        ("CCCCCCCCCCCCCCCCCC=C", "nonadec-1-ene"),
        ("CCCCCCCC=C", "1-nonene"),
        ("CCCC\\C=C/CC", "(Z)-3-octene"),
        ("[H]C(CC)=C([H])CCCC", "3-octene"),
        ("CCC=C", "but-1-ene"),
        ("C(CC)C(C)=C", "2-methyl-1-pentene"),
        ("CCCCCCC=C", "1-octene"),
        ("CCCC\\C=C\\CC", "(E)-3-octene"),
        ("CCCCCCCC\\C=C\\CCCCCCCC", "trans-octadec-9-ene"),
        ("[H]\\C(C)=C(/[H])C", "cis-but-2-ene"),
        ("CCC(C)=C", "2-methylbut-1-ene"),
        ("CCCCCCCCCCCCCCC=C", "1-hexadecene"),
        ("[H]C(CCCCCC)=C([H])CCCCCCCCCCCCCCC", "7-tricosene"),
        ("CC\\C=C/CC", "cis-3-hexene"),
        ("CCCCCCCCCC\\C=C/CCCCCC", "cis-octadec-7-ene"),
        ("CCCCCCCCC\\C=C/C", "(Z)-2-dodecene"),
        ("CCCCCCCCCCCCC=C", "1-pentadecene"),
        ("CCCCCCCCCCCCCCC=C", "1-heptadecene"),
        ("CCCCCCCCCCCCC\\C=C/CCCCCC", "(7Z)-tricosene"),
        ("CCCCCCCCCCC=C", "1-dodecene"),
        ("C(CCC)(CC(C)=C)C", "2,4-dimethyl-1-heptene"),
        ("C=C", "ethene"),
        ("CCCCCCCCCCCC\\C=C/CCCCCCCCC", "(10Z)-tricosene"),
        ("CC(C)C(C)CC=C", "4,5-dimethyl-1-hexene"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCC=C", "pentacos-1-ene"),
        ("CCCCCCCCC=C", "1-decene"),
        ("[H]C(CC)=C([H])CC", "3-hexene"),
        ("CC=C(C)C", "2-methylbut-2-ene"),
        ("CCCCCCCCCCCCCCCCCCCC=C", "1-docosene"),
        ("[H]C(CCCCCCCC)=C([H])CCCCCCCC", "octadec-9-ene"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCC=C", "nonacos-1-ene"),
        ("CC(C)CC(C)CC(C)CC=C", "4,6,8-trimethyl-1-nonene"),
        ("CCCCCC(C)CCC=C", "5-methyldec-1-ene"),
        ("CCCCCCCCC\\C=C\\C", "(E)-2-dodecene"),
        ("CCCCCCCCCCCCCCCC=C", "octadec-1-ene"),
        ("CCCCCCCCCC\\C=C\\C", "(2E)-2-tridecene"),
        ("C=CCCCCCC(CCC)C", "8-methyl-1-undecene"),
        ("CCCCCCCC\\C=C/CCCCCCCC", "cis-octadec-9-ene"),
        ("[H]C(C)=C([H])CCCCCCCCCCCCCCC", "octadec-2-ene"),
        ("CCCCCCCCCC=C", "1-undecene"),
        ("CCCC\\C=C/C", "(Z)-2-octene"),
        ("CC\\C=C\\CC", "trans-3-hexene"),
        ("C(CCCC(CCCC(C=C)=C)C)(CCCC(C)C)C", "neophytadiene"),
        ("CC(C)=C", "2-methylprop-1-ene"),
        ("[H]C(CC)=C([H])CCCCCCCC", "3-dodecene")
    ]

    for smi, name in test_examples:
        result, reason = is_alkene(smi)
        print(f"{name} (SMILES: {smi}) -> {result}. Reason: {reason}")
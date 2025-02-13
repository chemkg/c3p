"""
Classifies: CHEBI:47923 tripeptide
"""
#!/usr/bin/env python
"""
Classifies: Tripeptide – any oligopeptide consisting of three amino-acid residues connected by peptide linkages.
Heuristic improvements:
  • Instead of a fixed SMARTS pattern that may catch side-chain amide bonds,
    we iterate over bonds and count a peptide (backbone) bond when:
     - A bond connects a carbon and a nitrogen (in either order) with a single bond.
     - The carbon (carbonyl carbon) has at least one double bond to an oxygen.
     - The nitrogen (amide N) has at least one neighboring carbon (an α–carbon) apart from the carbonyl.
  • We also perform a rough check on the molecular weight.
  
A linear tripeptide should have 2 backbone peptide bonds and a cyclic one 3.
Note: This is a heuristic and may miss some edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as a peptide made of 3 amino-acid residues connected by peptide bonds.
    
    Heuristic applied:
      1. Count peptide bonds by iterating over bonds.
         A bond is considered a peptide bond if:
           - It connects a carbon (C) and a nitrogen (N) via a SINGLE bond.
           - The carbon (potential carbonyl) has a neighboring oxygen with a DOUBLE bond.
           - The nitrogen (amide N) has at least one neighbor (other than the carbonyl) that is a carbon.
      2. For a linear tripeptide, expect exactly 2 peptide bonds.
         For a cyclic tripeptide (where N- and C–termini are joined) expect exactly 3.
      3. Check that the molecular weight falls into a (very approximate) range.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a tripeptide, False otherwise.
        str: Explanation of the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count backbone peptide bonds by iterating over bonds.
    peptide_bond_count = 0
    for bond in mol.GetBonds():
        # We require a SINGLE bond between a carbon (atomic num 6) and nitrogen (atomic num 7).
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify which atom is carbon and which is nitrogen.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            carbon = a1
            nitrogen = a2
        elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 7:
            carbon = a2
            nitrogen = a1
        else:
            continue

        # Check that the carbon atom is a carbonyl carbon: it must have at least one double-bonded oxygen neighbor.
        has_carbonyl_oxygen = False
        for nbr in carbon.GetNeighbors():
            # Look for oxygen (atomic num 8) connected with a DOUBLE bond
            if nbr.GetAtomicNum() == 8:
                bond_to_nbr = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond_to_nbr is not None and bond_to_nbr.GetBondType() == Chem.BondType.DOUBLE:
                    has_carbonyl_oxygen = True
                    break
        if not has_carbonyl_oxygen:
            continue

        # Check the nitrogen atom: aside from the carbonyl carbon, it should be connected to at least one carbon
        # which we expect to be the α–carbon.
        has_alpha_carbon = False
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetIdx() == carbon.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                has_alpha_carbon = True
                break
        if not has_alpha_carbon:
            continue

        # If all conditions hold, count this bond as a backbone peptide bond.
        peptide_bond_count += 1

    # Check molecular weight – typical tripeptides are roughly 200-600 Da (allow some margin).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 1000):
        return False, f"Molecular weight {mol_wt:.1f} Da out of expected range for a typical tripeptide."

    # For a tripeptide we allow either:
    #  - Linear: exactly 2 peptide bonds (3 connected residues)
    #  - Cyclic: exactly 3 peptide bonds (the termini are linked)
    if peptide_bond_count == 2:
        return True, "Found 2 backbone peptide bonds consistent with a linear tripeptide (3 residues)."
    elif peptide_bond_count == 3:
        return True, "Found 3 backbone peptide bonds consistent with a cyclic tripeptide (3 residues)."
    else:
        return False, f"Found {peptide_bond_count} backbone peptide bond(s); expected 2 (linear) or 3 (cyclic) for a tripeptide."

# Example usage (testing a few SMILES strings):
if __name__ == "__main__":
    test_smiles = [
        # Glu-Phe-Val (linear): expected 2 backbone bonds.
        "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CCC(O)=O)C(O)=O",
        # Ala-Gly-His (linear): expected 2 backbone bonds.
        "C[C@H](N)C(=O)NCC(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"
    ]
    for smi in test_smiles:
        result, reason = is_tripeptide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
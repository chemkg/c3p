"""
Classifies: CHEBI:47923 tripeptide
"""
#!/usr/bin/env python
"""
Classifies: Tripeptide – any oligopeptide consisting of three amino-acid residues connected by peptide linkages.
This improved heuristic:
  • Iterates over bonds and counts only those backbone peptide bonds that satisfy stricter criteria:
      - A SINGLE bond connecting a carbon (atomic number 6) to a nitrogen (atomic number 7).
      - The carbon (carbonyl candidate) must have exactly one double-bonded oxygen neighbor and exactly three neighbors in total.
      - The nitrogen must have, apart from the carbonyl carbon, at least one carbon neighbor (as expected for an α–carbon).
  • From the candidate bonds the routine collects “backbone α–carbons” (neighbors of the peptide atoms, excluding the double bonded oxygen).
      In a genuine tripeptide the backbone will have exactly 3 such α–carbons.
  • Performs a rough molecular-weight check.

A linear tripeptide usually has 2 backbone peptide bonds whereas a cyclic one (where N- and C–termini join) has 3.
Note: This is a heuristic and may miss (or overcount) in heavily modified cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as 3 amino-acid residues connected by peptide bonds.
    
    Heuristic applied:
      1. Iterate over all bonds and identify candidate backbone peptide bonds.
         A bond is considered if:
           - It is a single bond connecting carbon (C) and nitrogen (N).
           - The carbon (potential carbonyl) has exactly one double-bonded oxygen neighbor
             and exactly three neighbors (one is the C=O, one is the peptide N, and one is the α–carbon).
           - The nitrogen (amide N) has at least one neighbor (aside from the carbonyl carbon) which is carbon.
      2. For each candidate bond, record the bond (by the carbon and nitrogen indices).
      3. Then, from these candidate bonds, collect the “backbone α–carbons” defined as any carbon atom
         that is directly attached to either the peptide carbonyl carbon (excluding the double-bonded oxygen)
         or to the amide nitrogen.
      4. A bona fide tripeptide should have exactly three such α–carbons.
      5. A further molecular weight filter (150–1000 Da) is applied.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tripeptide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    candidate_bonds = []  # list of tuples (carbon_idx, nitrogen_idx)

    # Loop over bonds to identify candidate backbone peptide bonds.
    for bond in mol.GetBonds():
        # Only consider SINGLE bonds.
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify the carbon and nitrogen atoms.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            carbon = a1
            nitrogen = a2
        elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 7:
            carbon = a2
            nitrogen = a1
        else:
            continue

        # Check that the carbon appears to be a carbonyl carbon:
        #   It must have at least one neighbor oxygen with a double bond.
        #   Also, we require exactly one such oxygen and that the carbon has exactly three neighbors.
        oxygen_double = False
        oxy_count = 0
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond_to_nbr = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond_to_nbr and bond_to_nbr.GetBondType() == Chem.BondType.DOUBLE:
                    oxygen_double = True
                    oxy_count += 1
        if not oxygen_double or oxy_count != 1:
            continue

        # Enforce that the carbonyl carbon has exactly three neighbors
        # (one is the double-bonded oxygen, one is the peptide nitrogen, and one is the α–carbon).
        if carbon.GetDegree() != 3:
            continue

        # Check the nitrogen: it must have at least one neighbor (other than the carbonyl carbon)
        # that is a carbon (expected to be the α–carbon).
        has_alpha_neighbor = False
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetIdx() == carbon.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                has_alpha_neighbor = True
                break
        if not has_alpha_neighbor:
            continue

        # All conditions met: count this bond as a candidate backbone peptide bond.
        candidate_bonds.append((carbon.GetIdx(), nitrogen.GetIdx()))
    
    # From the candidate bonds, collect the presumed backbone α–carbons.
    backbone_alpha = set()
    for (c_idx, n_idx) in candidate_bonds:
        carbon = mol.GetAtomWithIdx(c_idx)
        nitrogen = mol.GetAtomWithIdx(n_idx)
        # For the carbon (carbonyl): add any neighbor that is carbon and not the peptide nitrogen
        # and also not the oxygen (the one double-bonded to carbon).
        for nbr in carbon.GetNeighbors():
            if nbr.GetIdx() == n_idx:
                continue
            # Exclude the oxygen that is double-bonded.
            is_double_oxygen = False
            if nbr.GetAtomicNum() == 8:
                bond_to_nbr = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond_to_nbr and bond_to_nbr.GetBondType() == Chem.BondType.DOUBLE:
                    is_double_oxygen = True
            if not is_double_oxygen and nbr.GetAtomicNum() == 6:
                backbone_alpha.add(nbr.GetIdx())
        # For the nitrogen (amide N): add any neighbor that is carbon and not the carbonyl carbon.
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetIdx() == c_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                backbone_alpha.add(nbr.GetIdx())

    # A reasonable tripeptide (even with protecting groups) should have a backbone that gives exactly 3 α–carbons.
    if len(backbone_alpha) != 3:
        return False, f"Backbone α–carbon count is {len(backbone_alpha)}; expected 3 for a tripeptide."

    # Further filter by molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 1000):
        return False, f"Molecular weight {mol_wt:.1f} Da out of expected range for a typical tripeptide."

    # Finally, decide based on the number of candidate peptide bonds.
    # For a linear tripeptide we expect 2 backbone peptide bonds;
    # for a cyclic one (where the termini are linked) we expect 3.
    if len(candidate_bonds) == 2:
        return True, "Found 2 backbone peptide bonds consistent with a linear tripeptide (3 residues)."
    elif len(candidate_bonds) == 3:
        return True, "Found 3 backbone peptide bonds consistent with a cyclic tripeptide (3 residues)."
    else:
        # If extra candidate bonds appear (e.g. from protecting groups) but the backbone α–carbon count is correct,
        # we still classify the molecule as a tripeptide.
        return True, f"Backbone appears to have 3 residues with {len(candidate_bonds)} candidate peptide bonds."

# Example usage: testing a few SMILES strings.
if __name__ == "__main__":
    test_smiles = [
        # Glu-Phe-Val (expected linear tripeptide)
        "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CCC(O)=O)C(O)=O",
        # Ala-Gly-His (expected linear tripeptide)
        "C[C@H](N)C(=O)NCC(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"
    ]
    for smi in test_smiles:
        result, reason = is_tripeptide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
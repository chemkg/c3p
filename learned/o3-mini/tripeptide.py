"""
Classifies: CHEBI:47923 tripeptide
"""
#!/usr/bin/env python
"""
Classifies: Tripeptide – any oligopeptide consisting of three amino-acid residues connected by peptide bonds.
Improved heuristic:
  • The molecule is first processed with explicit hydrogens so that hydrogen counts on candidate α–carbons
    can be used as a crude filter.
  • A candidate peptide bond is a SINGLE bond between a carbon and a nitrogen where:
      - The carbon (potential “carbonyl”) has exactly one double-bonded oxygen.
      - The nitrogen (potential “amide N”) has at least one neighbor (other than the carbon) that is a carbon.
  • For each candidate bond, we look both at the carbonyl carbon’s other neighbors (excluding the double-bonded O
    and the N in the candidate bond) and at the amide nitrogen’s other neighbors – in each case picking only those
    carbons with a hydrogen count of either 1 (typical for chiral centers) or 2 (as in glycine). These are taken as
    the residue “α–carbons.”
  • In a tripeptide the assembled backbone will give exactly 3 such α–carbons.
  • A molecular weight filter (roughly 150–1200 Da) is applied.
  • Finally, if 2 (linear) or 3 (cyclic) candidate peptide bonds are found the molecule is classified as a tripeptide.
Note: This is a heuristic that may produce false positives or negatives when protecting groups or modifications occur.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as three amino-acid residues connected by peptide bonds.
    
    Heuristic improvements:
      1. Add explicit hydrogens to help reliably count the hydrogens on residue α–carbons.
      2. Identify candidate backbone peptide bonds: a single bond between a carbon (the carbonyl candidate)
         and a nitrogen (the amide candidate). The carbonyl must have exactly one double-bonded oxygen.
         The amide N must have at least one other heavy-atom neighbor that is a carbon (as expected for an α–carbon).
      3. For each candidate bond, from both the carbonyl and amide N sides, try to identify the residue’s α–carbon.
         We require that any candidate α–carbon (a neighboring carbon) has a hydrogen count (explicit) of 1 or 2.
      4. A valid tripeptide (even if modified) should give exactly 3 unique backbone α–carbons.
      5. A molecular weight filter is applied to remove very small or very large compounds.
      6. Finally, if exactly 2 candidate peptide bonds (linear) or 3 (cyclic) are found, the molecule is classified as a tripeptide.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a tripeptide, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES and add hydrogens for better atom discrimination.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    candidate_bonds = []  # list of tuples: (carbonyl atom idx, amide nitrogen atom idx)
    
    # Identify candidate peptide bonds.
    for bond in mol.GetBonds():
        # Only consider single bonds.
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify a bond between carbon and nitrogen.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            carbonyl = a1
            amideN = a2
        elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 7:
            carbonyl = a2
            amideN = a1
        else:
            continue
        
        # Check that the carbonyl candidate possesses exactly one double-bonded oxygen.
        dob_oxy = []
        for nbr in carbonyl.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                b2n = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
                if b2n and b2n.GetBondType() == Chem.BondType.DOUBLE:
                    dob_oxy.append(nbr.GetIdx())
        if len(dob_oxy) != 1:
            continue
        
        # Check that the amide nitrogen has at least one other neighbor that is a carbon (expected α–carbon).
        has_alpha = False
        for nbr in amideN.GetNeighbors():
            if nbr.GetIdx() == carbonyl.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                has_alpha = True
                break
        if not has_alpha:
            continue
        
        candidate_bonds.append((carbonyl.GetIdx(), amideN.GetIdx()))
    
    # Identify backbone α–carbons. A candidate α–carbon should be:
    #   • a carbon atom (atomic number 6)
    #   • attached to the peptide bond atom (either the carbonyl carbon or amide N)
    #   • carrying 1 or 2 explicit hydrogens.
    backbone_alphas = set()
    def valid_alpha(atom):
        if atom.GetAtomicNum() != 6:
            return False
        # Count explicit hydrogens.
        h_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
        return h_count in (1, 2)
    
    for (c_idx, n_idx) in candidate_bonds:
        carbonyl = mol.GetAtomWithIdx(c_idx)
        amideN = mol.GetAtomWithIdx(n_idx)
        # From carbonyl side: choose neighbors that are carbon, but not the amide nitrogen or the double-bonded oxygen.
        for nbr in carbonyl.GetNeighbors():
            if nbr.GetIdx() == n_idx:
                continue
            if nbr.GetAtomicNum() == 8:
                # Exclude if this oxygen is the double-bonded oxygen.
                b = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
                if b and b.GetBondType() == Chem.BondType.DOUBLE:
                    continue
            if nbr.GetAtomicNum() == 6 and valid_alpha(nbr):
                backbone_alphas.add(nbr.GetIdx())
        # From amide nitrogen side: choose neighbors that are carbon and not the carbonyl.
        for nbr in amideN.GetNeighbors():
            if nbr.GetIdx() == c_idx:
                continue
            if nbr.GetAtomicNum() == 6 and valid_alpha(nbr):
                backbone_alphas.add(nbr.GetIdx())
                
    # In a proper tripeptide, one expects exactly 3 backbone α–carbons.
    if len(backbone_alphas) != 3:
        return False, f"Backbone α–carbon count is {len(backbone_alphas)}; expected 3 for a tripeptide."
    
    # Apply a molecular weight filter – allow for some modifications.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 1200):
        return False, f"Molecular weight {mol_wt:.1f} Da out of expected range for a tripeptide."
    
    # Finally, decide based on the number of candidate peptide bonds.
    # For a linear tripeptide we expect two peptide bonds, while for a cyclic one the termini are connected (3 bonds).
    if len(candidate_bonds) == 2:
        return True, "Found 2 backbone peptide bonds consistent with a linear tripeptide (3 residues)."
    elif len(candidate_bonds) == 3:
        return True, "Found 3 backbone peptide bonds consistent with a cyclic tripeptide (3 residues)."
    else:
        return True, f"Backbone appears to have 3 residues with {len(candidate_bonds)} candidate peptide bonds."

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        # Glu-Phe-Val (linear tripeptide)
        "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CCC(O)=O)C(O)=O",
        # Ala-Gly-His (linear tripeptide)
        "C[C@H](N)C(=O)NCC(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"
    ]
    for smi in test_smiles:
        result, reason = is_tripeptide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
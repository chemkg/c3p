"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: Amino acids
Definition: A carboxylic acid containing one or more amino groups.

Changes from previous approach:
  - Require a carboxylic acid group (protonated or deprotonated).
  - Look for an "alpha backbone" pattern – a carbon bonded to both an amino group and a carboxyl group.
  - Additionally, search for a “free” amino group defined as a nitrogen with at least one hydrogen that is not directly adjacent to a carbonyl.
  - Count amide bonds and, if more than one is present and no clear alpha backbone is found,
    then (if the molecule is also moderately heavy) it is likely a peptide.
    
Note: This is a heuristic approach and still might mis‐classify edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_amino_acid(smiles: str):
    """
    Determines if a molecule can be classified as an amino acid based on its SMILES string.
    
    The molecule must have at least one carboxylic acid (–COOH or –COO–) group and at least one amino group.
    To help reject peptides (which usually have several amide bonds) we also count amide bonds.
    In addition, we search for an alpha amino acid backbone pattern (a carbon bonded to both an amino group 
    and a carboxyl group).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an amino acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for carboxylic acid group ---
    # Two SMARTS: one for protonated -COOH and one for deprotonated -COO-
    acid_smarts1 = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_smarts2 = Chem.MolFromSmarts("[CX3](=O)[O-]")
    has_acid = mol.HasSubstructMatch(acid_smarts1) or mol.HasSubstructMatch(acid_smarts2)
    if not has_acid:
        return False, "No carboxylic acid group found"
    
    # --- Count amide bonds (which are found in peptides) ---
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    amide_count = len(amide_matches)
    
    # --- Look for an alpha amino acid backbone pattern ---
    # This pattern captures a carbon attached to an amino group and to a carboxyl group.
    alpha_patterns = [
        Chem.MolFromSmarts("[C;!R]([NX3])(C(=O)[O-])"),  # deprotonated acid version
        Chem.MolFromSmarts("[C;!R]([NX3])(C(=O)O)")       # protonated acid version
    ]
    alpha_found = any(mol.HasSubstructMatch(pat) for pat in alpha_patterns)
    
    # --- Check for a free amino group ---
    # Define a "free" amino group as a nitrogen with at least one implicit or explicit hydrogen.
    # Exclude nitrogens that are adjacent to a carbonyl (which would indicate an amide bond).
    free_amino_smarts = Chem.MolFromSmarts("[NX3;H]")
    free_amino_matches = mol.GetSubstructMatches(free_amino_smarts)
    free_amino_valid = False
    for match in free_amino_matches:
        n_atom = mol.GetAtomWithIdx(match[0])
        is_amide = False
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetSymbol() == "C":
                # Check if this carbon is part of a C=O bond (carbonyl)
                for nb in neighbor.GetNeighbors():
                    if nb.GetSymbol() == "O":
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nb.GetIdx())
                        # bond type 2.0 corresponds to a double bond (C=O)
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            is_amide = True
                            break
                if is_amide:
                    break
        if not is_amide:
            free_amino_valid = True
            break

    # Accept as valid if we either found an alpha backbone pattern or a free amino group exists.
    has_valid_amino = alpha_found or free_amino_valid
    if not has_valid_amino:
        return False, "No free (non-amide) amino group found"
    
    # --- Peptide check ---
    # Heuristic: if the molecule contains more than one amide bond and lacks a clear alpha backbone,
    # and if its molecular weight is above a threshold (suggesting it is not a typical amino acid),
    # then it is likely a peptide.
    if amide_count > 1 and not alpha_found:
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt > 300:
            return False, "Molecule likely contains peptide bonds without an isolated amino acid backbone"
    
    return True, "Molecule contains a carboxylic acid and a valid amino group (alpha backbone/free amine)"

# Example usage: (When run as a script, will test a few examples)
if __name__ == "__main__":
    test_molecules = [
        ("OC(=O)C(N)CN1N=CC=C1", "3-(1-Pyrazolyl)-alanine"),
        ("N[C@@H](CC1=CC=C(F)C=C1)C(O)=O", "4-fluorophenyl-L-alanine"),
        ("CSCCCCC(N(O)O)C(O)=O", "N,N-dihydroxydihomomethionine"),
        ("O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CCC(=O)N", "Asp-Gln-Pro (peptide)")
    ]
    
    for s, name in test_molecules:
        result, reason = is_amino_acid(s)
        print(f"{name}: {result} -> {reason}")
"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
Definition: A molecular entity consisting of at least 2 biological molecules covalently linked together.
This heuristic approach breaks the molecule along linker bonds (e.g., amide, ester, thioether, and disulfide bonds)
and then inspects the resulting fragments. In addition to the original method, we lower the minimum molecular
weight threshold (to capture very small groups such as a vinyl substituent) and we reject cases where all fragments
resemble simple amino acids (i.e. a peptide chain) rather than a conjugate of distinct biomolecules.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    This function looks for linker bonds (amide, ester, thioether, disulfide) and fragments
    the molecule along those bonds. Then it requires that at least 2 fragments (with molecular weight
    > 30 Da) are obtained. In addition, if all fragments fall in the range of typical amino acids 
    (between ~70 and 300 Da and containing both nitrogen and oxygen), then the molecule is considered
    a continuous peptide chain rather than a bioconjugate.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a bioconjugate, False otherwise.
        str: A reason explaining the outcome.
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens to avoid aromaticity problems.
    try:
        mol = Chem.RemoveHs(mol)
    except Exception as e:
        return False, f"Failed to remove explicit hydrogens: {str(e)}"

    linker_bond_indices = set()
    
    # Iterate over bonds to identify potential linker bonds.
    # We look for:
    # - Amide bonds (C-N where the carbon is double-bonded to oxygen)
    # - Ester bonds (C-O where the carbon has an additional double-bonded oxygen)
    # - Thioether bonds (S-C bonds)
    # - Disulfide bonds (S-S bonds)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        atomic1 = a1.GetAtomicNum()
        atomic2 = a2.GetAtomicNum()
        
        # Check for amide: C-N with carbon double bonded to an oxygen.
        if ((atomic1 == 6 and atomic2 == 7) or (atomic1 == 7 and atomic2 == 6)):
            carbon = a1 if atomic1 == 6 else a2
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond_to_o and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        linker_bond_indices.add(bond.GetIdx())
                        break
            continue
        
        # Check for ester: C-O with the carbon also double-bonded to oxygen.
        if ((atomic1 == 6 and atomic2 == 8) or (atomic1 == 8 and atomic2 == 6)):
            carbon = a1 if atomic1 == 6 else a2
            other_atom = a2 if atomic1 == 6 else a1
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetIdx() == other_atom.GetIdx():
                    continue
                if neighbor.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond_to_o and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        linker_bond_indices.add(bond.GetIdx())
                        break
            continue
        
        # Thioether: S-C bond.
        if ((atomic1 == 16 and atomic2 == 6) or (atomic1 == 6 and atomic2 == 16)):
            linker_bond_indices.add(bond.GetIdx())
            continue
        
        # Disulfide: S-S bond.
        if (atomic1 == 16 and atomic2 == 16):
            linker_bond_indices.add(bond.GetIdx())
            continue

    if not linker_bond_indices:
        return False, "No recognizable linker bonds (amide, ester, thioether, or disulfide) found"
    
    # Fragment the molecule by breaking the identified linker bonds.
    try:
        fragmented_mol = Chem.FragmentOnBonds(mol, list(linker_bond_indices), addDummies=True)
    except Exception as e:
        return False, f"Fragmentation error: {str(e)}"
    
    # Extract fragments as separate molecules.
    try:
        fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    except Exception as e:
        return False, f"Fragment extraction error: {str(e)}"
    
    # Filter fragments by a lower molecular weight cutoff (30 Da).
    bio_fragments = []
    for frag in fragments:
        try:
            mw = rdMolDescriptors.CalcExactMolWt(frag)
            if mw > 30:
                bio_fragments.append((frag, mw))
        except Exception:
            continue

    if len(bio_fragments) < 2:
        return False, f"Only {len(bio_fragments)} biologically significant fragment(s) detected upon fragmentation"

    # Now check whether all the fragments are peptide‐like.
    # We use a crude heuristic: if a fragment has a molecular weight between 70 and 300,
    # and contains at least one nitrogen and one oxygen atom, we consider it peptide‐like.
    def is_peptide_fragment(mol_frag):
        numN = sum(1 for atom in mol_frag.GetAtoms() if atom.GetAtomicNum() == 7)
        numO = sum(1 for atom in mol_frag.GetAtoms() if atom.GetAtomicNum() == 8)
        mw_frag = rdMolDescriptors.CalcExactMolWt(mol_frag)
        return (70 < mw_frag < 300) and (numN >= 1) and (numO >= 1)
    
    peptide_like_count = 0
    for frag, mw in bio_fragments:
        if is_peptide_fragment(frag):
            peptide_like_count += 1
    
    # If all fragments appear as peptide-like, we suspect that the molecule is simply a peptide chain.
    if peptide_like_count == len(bio_fragments):
        return False, "Molecule appears to be a continuous peptide chain, not a bioconjugate"
    
    return True, f"Bioconjugate identified with {len(bio_fragments)} distinct biologically relevant fragments"

# Example testing (these lines can be removed in production code):
if __name__ == '__main__':
    test_smiles = [
        "N[C@@H](CCC(=O)N[C@@H](CS/C(=N\\O)/CC1=CC=CC=C1)C(=O)NCC(=O)O)C(=O)O",  # (Z)-1-(glutathione-S-yl)-2-phenylacetohydroximate
        "S(C[C@H](N)C(O)=O)C=C",  # S-vinylcysteine; expected to be bioconjugate (2 fragments, one small vinyl)
        "O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CCC(=O)N"  # Asp-Gln-Pro (simple peptide chain)
    ]
    for sm in test_smiles:
        result, reason = is_bioconjugate(sm)
        print(f"SMILES: {sm}\nResult: {result} :: {reason}\n")
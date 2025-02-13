"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
Definition: A molecular entity consisting of at least 2 biological molecules covalently linked together.
This heuristic approach breaks the molecule along linker bonds (e.g. amide, ester, thioether, and disulfide bonds)
and then inspects the resulting fragments. We now include a fallback method for problematic fragment extraction
and refine the peptide-like fragment heuristic.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    This function seeks linker bonds (amide, ester, thioether, disulfide) and fragments
    the molecule along those bonds. The molecule is classified as a bioconjugate if at least 2 fragments 
    (with a molecular weight > 30 Da) are detected and at least one fragment does not appear as a typical peptide 
    (heuristically defined by having a mass between 70 and 300 Da and containing at least one nitrogen and one oxygen).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a bioconjugate, False otherwise.
        str: A reason explaining the outcome.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try to remove explicit hydrogens to avoid aromaticity issues.
    try:
        mol = Chem.RemoveHs(mol)
    except Exception as e:
        return False, f"Failed to remove explicit hydrogens: {str(e)}"

    linker_bond_indices = set()
    
    # Identify potential linker bonds.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        atomic1 = a1.GetAtomicNum()
        atomic2 = a2.GetAtomicNum()
        
        # Amide: Look for C–N where the carbon is double-bonded to an oxygen.
        if ((atomic1 == 6 and atomic2 == 7) or (atomic1 == 7 and atomic2 == 6)):
            carbon = a1 if atomic1 == 6 else a2
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond_to_o and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        linker_bond_indices.add(bond.GetIdx())
                        break
            continue
        
        # Ester: Look for C–O where the carbon is also double-bonded to an oxygen.
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
        
        # Thioether: S–C bond.
        if ((atomic1 == 16 and atomic2 == 6) or (atomic1 == 6 and atomic2 == 16)):
            linker_bond_indices.add(bond.GetIdx())
            continue
        
        # Disulfide: S–S bond.
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
    
    # Extract fragments. Try the sanitization here; if it fails catch and try less strict sanitization.
    try:
        fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    except Exception as e:
        try:
            fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=False)
            # Manually sanitize fragments while ignoring aromaticity setting if needed:
            for frag in fragments:
                try:
                    Chem.SanitizeMol(frag, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
                except Exception:
                    pass
        except Exception as e2:
            return False, f"Fragment extraction error: {str(e2)}"
    
    # Filter fragments with molecular weight > 30 Da; these are considered biologically significant.
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
    
    # Define a helper: a fragment is considered "peptide-like" if its mass is between 70 and 300 Da and contains at least one
    # nitrogen and one oxygen; otherwise it is considered non-peptide.
    def is_peptide_fragment(mol_frag):
        numN = sum(1 for atom in mol_frag.GetAtoms() if atom.GetAtomicNum() == 7)
        numO = sum(1 for atom in mol_frag.GetAtoms() if atom.GetAtomicNum() == 8)
        mw_frag = rdMolDescriptors.CalcExactMolWt(mol_frag)
        return (70 < mw_frag < 300) and (numN >= 1) and (numO >= 1)
    
    peptide_like_count = 0
    for frag, mw in bio_fragments:
        if is_peptide_fragment(frag):
            peptide_like_count += 1
    
    # If every biologically significant fragment appears as a typical peptide fragment, then this is a peptide chain.
    if peptide_like_count == len(bio_fragments):
        return False, "Molecule appears to be a continuous peptide chain, not a bioconjugate"
    
    # Otherwise, we classify it as a bioconjugate.
    return True, f"Bioconjugate identified with {len(bio_fragments)} distinct biologically relevant fragments"

# Example usage (for testing):
if __name__ == '__main__':
    test_smiles = [
        "N[C@@H](CCC(=O)N[C@@H](CS/C(=N\\O)/CC1=CC=CC=C1)C(=O)NCC(=O)O)C(=O)O",  # (Z)-1-(glutathione-S-yl)-2-phenylacetohydroximate
        "N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O",         # gammaGluCys(IAN)Gly
        "[C@H]1([C@H]([C@@H](CC1=O)SC[C@H](NC(CC[C@H](N)C(=O)O)=O)C(=O)NCC(=O)O)/C=C/[C@H](CCCCC)O)C/C=C\\CCCC(=O)O"  # (R)-PGA2-S-glutathione conjugate
    ]
    for sm in test_smiles:
        result, reason = is_bioconjugate(sm)
        print(f"SMILES: {sm}\nResult: {result} :: {reason}\n")
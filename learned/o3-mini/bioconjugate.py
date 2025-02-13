"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
Definition: A molecular entity consisting of at least 2 biological molecules covalently linked together.
This heuristic approach breaks the molecule along bonds connecting bio‐moieties (e.g., amide, ester, thioether, and disulfide bonds)
to see if it can be split into at least two biologically significant fragments.
Improved error handling is included to catch aromaticity or fragmentation issues.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    The method detects linker bonds (amide, ester, thioether, disulfide) and then fragments
    the molecule along these bonds. If at least 2 fragments with molecular weight > 50 Da are obtained,
    it is considered a bioconjugate.

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
    
    # Remove explicit hydrogens to avoid possible aromaticity issues.
    try:
        mol = Chem.RemoveHs(mol)
    except Exception as e:
        return False, f"Failed to remove explicit hydrogens: {str(e)}"
    
    linker_bond_indices = []
    
    # Iterate over bonds to identify potential linker bonds:
    # Check for: amide bonds, ester bonds, thioether bonds, and disulfide bonds.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify amide bonds: a C-N bond where the carbon is doubly bonded to oxygen.
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7) or 
            (a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6)):
            carbon = a1 if a1.GetAtomicNum() == 6 else a2
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond_to_o and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        linker_bond_indices.append(bond.GetIdx())
                        break
            continue
        
        # Identify ester bonds: a C-O bond where the carbon has an additional oxygen via a double bond.
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8) or 
            (a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6)):
            carbon = a1 if a1.GetAtomicNum() == 6 else a2
            other_atom = a2 if a1.GetAtomicNum() == 6 else a1
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetIdx() == other_atom.GetIdx():
                    continue
                if neighbor.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond_to_o and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        linker_bond_indices.append(bond.GetIdx())
                        break
            continue
        
        # Identify thioether bonds: S-C bond.
        if ((a1.GetAtomicNum() == 16 and a2.GetAtomicNum() == 6) or 
            (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 16)):
            linker_bond_indices.append(bond.GetIdx())
            continue
        
        # Identify disulfide bonds: S-S bond.
        if (a1.GetAtomicNum() == 16 and a2.GetAtomicNum() == 16):
            linker_bond_indices.append(bond.GetIdx())
            continue

    linker_bond_indices = list(set(linker_bond_indices))
    if not linker_bond_indices:
        return False, "No recognizable linker bonds (amide, ester, thioether, or disulfide) found"

    # Fragment the molecule by breaking identified bonds.
    try:
        fragmented_mol = Chem.FragmentOnBonds(mol, linker_bond_indices, addDummies=True)
    except Exception as e:
        return False, f"Fragmentation error: {str(e)}"
    
    # Obtain the fragments as separate molecules.
    try:
        # sanitizeFrags=True attempts to clean each fragment.
        fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    except Exception as e:
        return False, f"Fragment extraction error: {str(e)}"
    
    # Filter fragments that are of biological relevance (molecular weight > 50 Da).
    bio_fragments = []
    for frag in fragments:
        try:
            mw = rdMolDescriptors.CalcExactMolWt(frag)
            if mw > 50:
                bio_fragments.append(frag)
        except Exception:
            continue

    if len(bio_fragments) < 2:
        return False, f"Only {len(bio_fragments)} biologically significant fragment(s) detected upon fragmentation"
    
    return True, f"Bioconjugate identified with {len(bio_fragments)} distinct biologically relevant fragments"
"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
Definition: A molecular entity consisting of at least 2 biological molecules covalently linked together.
This heuristic approach fragments the molecule along bonds often connecting bio‐moieties.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    Our approach looks for “linker” bonds (such as amide, ester, thioether and disulfide bonds)
    whose cleavage yields at least two fragments of significant size.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    linker_bond_indices = []
    
    # Iterate over bonds to identify potential bioconjugation bonds:
    # Check for: amide bonds, ester bonds, thioether bonds, and disulfide bonds.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify amide bonds: look for a C-N bond in which the carbon has a carbonyl (double bonded to O)
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7) or (a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6)):
            # Identify which atom is carbon
            carbon = a1 if a1.GetAtomicNum() == 6 else a2
            # Check neighbors of the carbon for a double-bonded oxygen (the carbonyl)
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond_to_o and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        linker_bond_indices.append(bond.GetIdx())
                        break
            continue  # continue with next bond if amide identified
        
        # Identify ester bonds: a C-O bond with the carbon bound to an additional oxygen via a double bond 
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8) or (a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6)):
            carbon = a1 if a1.GetAtomicNum() == 6 else a2
            other_atom = a2 if a1.GetAtomicNum() == 6 else a1
            # Look for a carbonyl O attached to the same carbon (other than the current oxygen)
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetIdx() == other_atom.GetIdx():
                    continue
                if neighbor.GetAtomicNum() == 8:
                    bond_to_o = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond_to_o and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
                        linker_bond_indices.append(bond.GetIdx())
                        break
            continue
        
        # Identify thioether bonds: S-C bond
        if ((a1.GetAtomicNum() == 16 and a2.GetAtomicNum() == 6) or (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 16)):
            linker_bond_indices.append(bond.GetIdx())
            continue
        
        # Identify disulfide bonds: S-S bond
        if (a1.GetAtomicNum() == 16 and a2.GetAtomicNum() == 16):
            linker_bond_indices.append(bond.GetIdx())
            continue

    # Remove duplicate bond indices if any
    linker_bond_indices = list(set(linker_bond_indices))
    
    if not linker_bond_indices:
        return False, "No recognizable linker bonds (amide, ester, thioether, or disulfide) found"
    
    # Fragment the molecule by "breaking" the identified bonds.
    # The addDummies flag will add dummy atoms at the break points.
    fragmented_mol = Chem.FragmentOnBonds(mol, linker_bond_indices, addDummies=True)
    # Get the fragments as separate molecules
    fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    
    # For each fragment, we can compute its molecular weight.
    # We consider only fragments with a molecular weight > 50 Da as biologically significant.
    bio_fragments = []
    for frag in fragments:
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        if mw > 50:
            bio_fragments.append(frag)
            
    # A bioconjugate is defined as having at least 2 biologically significant fragments.
    if len(bio_fragments) < 2:
        return False, f"Only {len(bio_fragments)} biologically significant fragment(s) detected upon fragmentation"
    
    return True, f"Bioconjugate identified with {len(bio_fragments)} distinct biologically relevant fragments"

# For testing, you might run (uncomment the following lines):
# test_smiles = "N[C@@H](CCC(=O)N[C@@H](CS/C(=N\\O)/CC1=CC=CC=C1)C(=O)NCC(=O)O)C(=O)O"  # (Z)-1-(glutathione-S-yl)-2-phenylacetohydroximate
# result, reason = is_bioconjugate(test_smiles)
# print(result, reason)
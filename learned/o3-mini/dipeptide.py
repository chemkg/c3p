"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: Any molecule that contains two amino‐acid residues connected by peptide linkages.
The algorithm here first searches for two candidate α–carbons (the central carbon in an amino acid residue) 
using a heuristic: an sp3 carbon that is directly bonded to at least one nitrogen (usually from an amide bond) 
and directly bonded to a carbon that is “carbonyl‐like” (i.e. that has a C=O double bond). Then the code 
checks that one of those carbonyl carbons is also bonded to a nitrogen that in turn is bonded (on its other side)
to the other α–carbon. This connection is taken as evidence for a peptide bond linking the two amino acid residues.
Note: Because dipeptides are often modified (e.g. with acyl groups on the N‐terminus or cyclized backbones), 
this method is heuristic and may not work in all cases. If too many ambiguities are present the function may be revised.
"""

from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is defined as any molecule that contains two amino‐acid residues connected by peptide linkage(s).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a dipeptide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic: Identify candidate alpha carbons.
    # We define an alpha carbon as any (implicit‐H) carbon atom (atomic number 6) that is bonded to at least one nitrogen
    # and also bonded to at least one carbon that is engaged in a carbonyl (C=O) bond.
    alpha_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        neighbors = atom.GetNeighbors()
        # Check for at least one nitrogen neighbor.
        has_nitrogen = any(nb.GetAtomicNum() == 7 for nb in neighbors)
        # Check for at least one neighboring carbon that is carbonyl-like:
        has_carbonyl = False
        for nb in neighbors:
            if nb.GetAtomicNum() == 6:
                # Look at bonds of neighbor to see if there is a double bond to oxygen.
                for bond in nb.GetBonds():
                    # Ensure the bond involves an oxygen (atomic number 8) and is a double bond.
                    if bond.GetBondTypeAsDouble() == 2.0:
                        begin = bond.GetBeginAtom()
                        end = bond.GetEndAtom()
                        if (begin.GetAtomicNum() == 8 or end.GetAtomicNum() == 8):
                            has_carbonyl = True
                            break
                if has_carbonyl:
                    break
        if has_nitrogen and has_carbonyl:
            alpha_indices.append(atom.GetIdx())
    
    if len(alpha_indices) != 2:
        return False, f"Found {len(alpha_indices)} candidate alpha–carbons; expected exactly 2 for a dipeptide"
    
    # Now check that these two candidate alpha–carbons are connected via a peptide linkage.
    # For a typical dipeptide, one amino acid residue (let’s call it residue 1) has its alpha carbon connected to its own carboxyl carbon.
    # That carboxyl carbon in turn is bonded (via an amide bond) to a nitrogen that is attached to the second residue’s alpha carbon.
    found_peptide_link = False
    for a1_idx in alpha_indices:
        a1 = mol.GetAtomWithIdx(a1_idx)
        # Look over neighbors of a1 to find a carbon that could be the carbonyl in residue 1.
        for nb in a1.GetNeighbors():
            if nb.GetAtomicNum() == 6:
                # Test if nb has a double bond to oxygen (a carbonyl group)
                carbonyl_found = False
                for bond in nb.GetBonds():
                    # bond between nb and oxygen with double bond
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nb)
                        if other.GetAtomicNum() == 8:
                            carbonyl_found = True
                            break
                if not carbonyl_found:
                    continue
                # Now nb is a carbonyl carbon candidate that is bonded to a1.
                # Look at its neighbors (other than a1) for an amide nitrogen.
                for nb2 in nb.GetNeighbors():
                    if nb2.GetIdx() == a1_idx:
                        continue
                    if nb2.GetAtomicNum() == 7:  # candidate amide nitrogen
                        # Now check if nb2 is bonded to the other alpha carbon.
                        for nb3 in nb2.GetNeighbors():
                            if nb3.GetIdx() == nb.GetIdx():
                                continue
                            if nb3.GetIdx() in alpha_indices and nb3.GetIdx() != a1_idx:
                                found_peptide_link = True
                                break
                        if found_peptide_link:
                            break
                if found_peptide_link:
                    break
        if found_peptide_link:
            break

    if not found_peptide_link:
        return False, "No peptide linkage connecting two amino–acid residues found"
    
    return True, "Found two amino–acid residues connected by a peptide linkage"

# Example usage (uncomment the following lines to test):
# test_smiles = "SC[C@H](NC(=O)[C@@H](N)CCCCN)C(O)=O"  # Lys-Cys example
# result, reason = is_dipeptide(test_smiles)
# print(result, reason)
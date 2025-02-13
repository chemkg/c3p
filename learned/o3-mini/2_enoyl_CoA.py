"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
Definition: An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond 
            between positions 2 and 3. Here position 1 is the thioester carbonyl carbon.
            
The algorithm first verifies that a CoA moiety is present (using a characteristic SMARTS fragment)
and then looks for a thioester carbonyl (a C=O that is directly bound to a sulfur atom).
Among those, it identifies the acyl (R) group attached to the carbonyl carbon (excluding the S),
and then checks whether the alpha carbon (the one directly attached) is involved in a double bond 
to a beta carbon (and that double bond does not go back to the carbonyl).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) qualifies as a 2-enoyl-CoA.
    To qualify the molecule must contain:
      1. A CoA moiety as detected by a CoA substructure fragment.
      2. An acyl thioester group in which the acyl chain (the group attached to the carbonyl C)
         has a double bond between its first (alpha) and second (beta) carbon.
         
    Args:
        smiles (str): SMILES string representation of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a 2-enoyl-CoA; False otherwise.
        str: An explanation for the decision.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Step 1. Check for the presence of a CoA moiety.
    # We use a characteristic fragment from the pantetheine/ADP portion of CoA.
    # (This SMARTS works reasonably well for many acyl-CoA structures.)
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in SMARTS for CoA fragment."
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety fragment not detected."

    # Step 2. Identify thioester carbonyls.
    # Look for carbon atoms that have:
    #   - a double bond to an oxygen (carbonyl)
    #   - a single bond to a sulfur (thioester linkage)
    # For each such carbon, check its other substituent (the acyl chain) for the expected unsaturation.
    enoyl_found = False

    for atom in mol.GetAtoms():
        # Only consider carbon atoms
        if atom.GetAtomicNum() != 6:
            continue

        # Look for a carbonyl oxygen (C=O) and a sulfur neighbor.
        has_carbonyl_oxygen = False
        has_S_neighbor = False
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_carbonyl_oxygen = True
            if nbr.GetAtomicNum() == 16 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                has_S_neighbor = True
        # Skip atoms that are not part of a thioester carbonyl.
        if not (has_carbonyl_oxygen and has_S_neighbor):
            continue

        # For a carbonyl carbon in a thioester, one neighbor is the sulfur (leading to CoA).
        # The other neighbor (if any) is the acyl chain (R group); that should be the alpha carbon.
        acyl_candidates = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() not in (8, 16):  # Exclude oxygen and sulfur
                # Also, ensure that the bond between the carbonyl and this atom is a single bond.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    acyl_candidates.append(nbr)
        # We expect a thioester carbonyl to have exactly one acyl substituent.
        if not acyl_candidates:
            continue
        # For each candidate (usually there is just one), check for a double bond on the candidate.
        for alpha in acyl_candidates:
            # The expected pattern: the alpha carbon (attached to the carbonyl) should be 
            # involved in a double bond with a beta carbon. To avoid picking up distant unsaturation,
            # we only consider double bonds emanating from the alpha carbon that do not go back to
            # the carbonyl atom.
            for bond in alpha.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    beta = bond.GetOtherAtom(alpha)
                    if beta.GetIdx() == atom.GetIdx():
                        continue  # Skip if it goes back to the carbonyl carbon.
                    # Found a double bond from the alpha carbon.
                    enoyl_found = True
                    break
            if enoyl_found:
                break
        if enoyl_found:
            break

    if not enoyl_found:
        return False, "No thioester acyl chain with a double bond between its alpha and beta carbons was detected."

    return True, "Molecule contains a 2-enoyl-CoA fragment (unsaturated acyl thioester with CoA moiety)."

# Example usage:
# test_smiles = "CCCCCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
# result, explanation = is_2_enoyl_CoA(test_smiles)
# print(result, explanation)
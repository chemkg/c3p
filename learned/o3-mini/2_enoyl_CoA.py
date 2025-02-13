"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
Definition: An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond 
            between positions 2 and 3. Here position 1 is the thioester carbonyl carbon.
            
This algorithm first verifies that the CoA moiety is present (using a characteristic SMARTS fragment).
Then, it searches for a thioester group – a carbonyl carbon double bonded to oxygen and singly bonded to a sulfur.
Using the CoA match we require that the sulfur is part of the known CoA fragment while the acyl (R) group is external.
Finally, among the acyl side–chain atoms (the candidate alpha carbon, i.e. the one attached to the carbonyl,
other than the S) we check for a double bond to a beta carbon (which does not go back to the carbonyl).
This is taken as evidence that the molecule contains a 2-enoyl-CoA fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) qualifies as a 2-enoyl-CoA.
    To qualify the molecule must contain:
      1. A CoA moiety (detected via a characteristic SMARTS fragment).
      2. A thioester carbonyl (C=O bonded to sulfur) in which the S is part of the CoA.
      3. An acyl chain (the R group attached to the carbonyl carbon that is NOT part of CoA)
         in which the alpha carbon (directly bonded to the carbonyl) carries a double bond 
         (to a beta carbon that is not the carbonyl).
         
    Args:
        smiles (str): SMILES string representation of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a 2-enoyl-CoA; False otherwise.
        str: An explanation for the decision.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1. Check for presence of CoA moiety.
    # We use a SMARTS fragment that covers a characteristic portion of the CoA (pantetheine/ADP part).
    # Note: This fragment is heuristic and works reasonably well for many acyl-CoA molecules.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in SMARTS for CoA fragment."
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "CoA moiety fragment not detected."
    
    # Form a set of all atom indices in any CoA match (to know which atoms belong to CoA).
    coa_atom_set = set()
    for match in coa_matches:
        coa_atom_set.update(match)
    
    # Step 2. Identify thioester carbonyl groups.
    # We loop over carbon atoms looking for a C=O (double bonded to O) and a single bond to S.
    # For each such carbon, we then pick the substituent that is not oxygen or the sulfur.
    # We further require that the sulfur neighbor is part of CoA (from our match)
    # while the acyl substituent is outside the CoA fragment.
    enoyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Look for the required bonds: a double bond to oxygen and a single bond to sulfur.
        has_carbonyl_oxygen = False
        s_neighbors = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_carbonyl_oxygen = True
            if nbr.GetAtomicNum() == 16 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                s_neighbors.append(nbr)
        if not (has_carbonyl_oxygen and s_neighbors):
            continue
        
        # For each sulfur neighbor, require that it is part of the known CoA fragment.
        coa_s = None
        for s in s_neighbors:
            if s.GetIdx() in coa_atom_set:
                coa_s = s
                break
        if coa_s is None:
            # This thioester carbonyl does not link to CoA; skip.
            continue
        
        # Now, among the other substituents of the carbonyl carbon, select the acyl chain candidate.
        acyl_candidates = []
        for nbr in atom.GetNeighbors():
            # Exclude the oxygen (carbonyl) and the sulfur (which is in CoA).
            if nbr.GetAtomicNum() in (8, 16):
                continue
            # Also require that the bond is a single bond.
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            # We also require that this candidate is not part of the CoA fragment.
            if nbr.GetIdx() in coa_atom_set:
                continue
            acyl_candidates.append(nbr)
            
        if not acyl_candidates:
            continue
        
        # For each candidate alpha carbon, check if it carries the required double bond.
        for alpha in acyl_candidates:
            # We look for at least one double bond from the alpha carbon that goes to a beta carbon.
            # Ensure that the beta is not the carbonyl carbon.
            for bond in alpha.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    beta = bond.GetOtherAtom(alpha)
                    if beta.GetIdx() == atom.GetIdx():
                        continue  # Skip the carbonyl back-bond.
                    # Found a double bond starting at the alpha carbon.
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
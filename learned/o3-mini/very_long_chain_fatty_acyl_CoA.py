"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: very long-chain fatty acyl-CoA
Definition: a fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.
In this implementation the acyl chain is extracted by identifying the thioester bond
(i.e. a carbonyl C with a double bond to an oxygen that is attached to a sulfur)
and then fragmenting the molecule at that bond. The fragment that does not contain sulfur
is assumed to be the fatty acyl portion. Its carbon count is then compared to 22.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a very long-chain fatty acyl-CoA.
    For our purposes, a molecule qualifies if it contains a thioester bond (C(=O)-S) linking an acyl group 
    whose number of carbon atoms is greater than 22.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the fatty acyl branch has more than 22 carbon atoms, False otherwise.
        str: Reason for the classification.
    """
    # Attempt to create a molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, search for a thioester bond: a bond connecting a carbon and a sulfur,
    # where the carbon is part of a carbonyl group (i.e. has a double bond to an oxygen).
    thioester_bond_idx = None
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # One atom must be sulfur (atomic num 16) and the other a carbon (atomic num 6)
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 16) or 
            (a1.GetAtomicNum() == 16 and a2.GetAtomicNum() == 6)):
            # Identify the carbon atom from the pair.
            carbon = a1 if a1.GetAtomicNum() == 6 else a2
            sulfur = a2 if a1.GetAtomicNum() == 6 else a1
            
            # Check if the carbon has a double bond to an oxygen.
            found_carbonyl = False
            for nb in carbon.GetBonds():
                # Exclude the bond to the sulfur.
                if nb.GetOtherAtom(carbon).GetIdx() == sulfur.GetIdx():
                    continue
                # Check if neighboring atom is oxygen with bond order 2.
                if nb.GetOtherAtom(carbon).GetAtomicNum() == 8 and nb.GetBondTypeAsDouble() >= 2.0:
                    found_carbonyl = True
                    break
            if found_carbonyl:
                thioester_bond_idx = bond.GetIdx()
                break

    if thioester_bond_idx is None:
        return False, "No thioester (C(=O)-S) bond found; not a fatty acyl-CoA structure"
    
    # Now break the molecule at the thioester bond to separate the acyl (fatty) fragment 
    # from the CoA portion. We'll use FragmentOnBonds with addDummies=False.
    fragMol = rdmolops.FragmentOnBonds(mol, [thioester_bond_idx], addDummies=False)
    
    # Get disconnected fragments
    frags = Chem.GetMolFrags(fragMol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "Unable to fragment the molecule"

    # The acyl chain is expected to be the fragment that does NOT contain sulfur.
    candidate_frags = []
    for frag in frags:
        # If the fragment has any sulfur atom, we assume it is part of the CoA moiety.
        if not any(atom.GetAtomicNum() == 16 for atom in frag.GetAtoms()):
            candidate_frags.append(frag)
    
    if not candidate_frags:
        return False, "No candidate acyl chain fragment found after breaking the thioester bond"
    
    # If there is more than one candidate, choose the fragment with the highest number of carbon atoms.
    acyl_frag = max(candidate_frags, key=lambda m: sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6))
    carbon_count = sum(1 for atom in acyl_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # For debugging, one might print the candidate fragment SMILES:
    # print("Acyl fragment SMILES:", Chem.MolToSmiles(acyl_frag))
    
    if carbon_count > 22:
        return True, f"Acyl chain fragment has {carbon_count} carbon atoms > 22"
    else:
        return False, f"Acyl chain fragment has only {carbon_count} carbon atoms (needs > 22)"

# Example use:
# result, reason = is_very_long_chain_fatty_acyl_CoA("CCCCCCCCC\\C=C/C(=O)SCCNC(=O)...")
# print(result, reason)
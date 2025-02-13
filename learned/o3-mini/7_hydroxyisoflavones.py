"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python3
"""
Classifies: 7-hydroxyisoflavones
A 7-hydroxyisoflavone is defined as a hydroxyisoflavone compound having a hydroxy (–OH)
group at the 7-position of the isoflavone A ring.
The algorithm uses a two‐step approach:
  1. Check that the molecule contains an isoflavone core. Here we use a SMARTS for a 3-phenylchromen-4-one.
  2. Check that at least one aromatic carbon in a six‐membered (benzene) ring fused to the pyranone
     bears a direct –OH substituent.
If either check fails, the molecule is not classified as a 7-hydroxyisoflavone.
Note: This procedure is heuristic and may not correctly classify all edge cases.
"""

from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    
    The method first verifies that the molecule has an isoflavone core (3-phenylchromen-4-one)
    and then looks for an –OH group on one of the aromatic rings of the fused bicyclic core,
    as expected for a hydroxyl at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS for an isoflavone core.
    # The pattern here is a simplified representation of 3-phenylchromen-4-one.
    # It requires a benzopyran-4-one fused system with a phenyl substituent.
    iso_core_smarts = "c1cc2oc(=O)c(c2c1)-c1ccccc1"
    iso_core = Chem.MolFromSmarts(iso_core_smarts)
    if iso_core is None:
        return False, "Error in SMARTS pattern."
        
    # Check if the isoflavone core is present.
    if not mol.HasSubstructMatch(iso_core):
        return False, "Molecule does not contain the expected isoflavone core."
    
    # Get ring information for further analysis
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Look for hydroxyl groups (-OH) on aromatic carbons in a benzene ring that is fused in the core.
    # We iterate over atoms: if an aromatic carbon has a neighboring oxygen that is an –OH, 
    # and if that carbon is in a six-membered ring that also contains a carbonyl group (C=O),
    # we assume that this –OH is on the A ring (typical of a 7-hydroxyisoflavone).
    for atom in mol.GetAtoms():
        # Check for aromatic carbon
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            # Look at neighbors to see if one is an -OH group (oxygen atom with at least one H)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    # Now check in which rings this atom occurs.
                    for ring in ring_info:
                        if atom.GetIdx() in ring and len(ring) == 6:  # Likely a benzene ring
                            # Look in the ring for a carbonyl.
                            has_carbonyl = False
                            for idx in ring:
                                ring_atom = mol.GetAtomWithIdx(idx)
                                # Examine neighbors of ring_atom to see if one is a double-bonded O (C=O)
                                for neighbor in ring_atom.GetNeighbors():
                                    bond = mol.GetBondBetweenAtoms(ring_atom.GetIdx(), neighbor.GetIdx())
                                    # bond type check: double bond and oxygen
                                    if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                                        has_carbonyl = True
                                        break
                                if has_carbonyl:
                                    break
                            if has_carbonyl:
                                # Found an aromatic carbon (in a benzene ring of the fused core)
                                # that bears an -OH and where the ring also contains a carbonyl.
                                return True, "Molecule has isoflavone core with a hydroxyl on the A ring (7-position likely)."
    
    return False, "Molecule has isoflavone core but no hydroxyl found on the A ring (7-position)."


# For testing the function with one example (uncomment below lines to run a quick test)
# if __name__ == "__main__":
#     test_smiles = "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O"  # 7-hydroxyisoflavone
#     result, reason = is_7_hydroxyisoflavones(test_smiles)
#     print(result, reason)
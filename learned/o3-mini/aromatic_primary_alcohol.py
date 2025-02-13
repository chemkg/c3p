"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Primary Alcohol
Definition: Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon 
            which is itself bonded to an aromatic ring.
            
Improvement notes:
 – Instead of iterating over OH atoms, we loop over carbons. For each carbon,
   we check that it is sp3 and carries exactly two hydrogens (i.e. CH2).
 – We then check that it is singly bonded to exactly one oxygen that qualifies as a hydroxyl group,
   and that its one remaining heavy-atom neighbor is aromatic.
 – This minimizes false positives (where extra substituents might be present) and avoids missing groups
   because of implicit hydrogens/stereochemistry.
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains an aromatic primary alcohol group.
    A qualifying structure contains a CH2OH group (primary alcohol) 
    where the CH2 carbon is sp3, bears exactly two hydrogens, has a single exocyclic -OH (via a single bond),
    and its other (and only) heavy-atom neighbor is aromatic.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule contains an aromatic primary alcohol group, False otherwise.
        str: Explanation/reason for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are accurate
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms looking for a candidate carbon for primary alcohol
    for atom in mol.GetAtoms():
        # We want a carbon
        if atom.GetAtomicNum() != 6:
            continue
        # Ensure the carbon is sp3 (primary alcohol carbons are sp3)
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # Check it has exactly two hydrogens (which implies a CH2 group)
        if atom.GetTotalNumHs() != 2:
            continue
        
        # Identify oxygen neighbors connected via single bonds that might be the -OH group.
        hydroxyl_oxygens = []
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check the bond between the atom and oxygen is a single bond (avoid carbonyl O, etc.)
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # Check that the oxygen is indeed -OH:
                # It should have at least one hydrogen neighbor.
                h_neighbors = [n for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 1]
                if len(h_neighbors) < 1:
                    continue
                hydroxyl_oxygens.append(neighbor)
                
        # We require exactly one hydroxyl oxygen connected to the carbon
        if len(hydroxyl_oxygens) != 1:
            continue
        
        # Now, the carbon should have exactly one other heavy-atom neighbor apart from the hydroxyl oxygen.
        # (Primary alcohol: R-CH2OH where R is one substituent)
        other_neighbors = [n for n in atom.GetNeighbors() 
                           if n.GetAtomicNum() > 1 and n.GetIdx() != hydroxyl_oxygens[0].GetIdx()]
        if len(other_neighbors) != 1:
            continue
        
        # Check that this single substituent is an aromatic atom.
        # (This ensures that the –CH2– is directly attached to an aromatic ring.)
        if not other_neighbors[0].GetIsAromatic():
            continue
        
        # If all conditions are met, we have found an aromatic primary alcohol group.
        return True, "Found an aromatic primary alcohol group (CH2OH attached to an aromatic ring)"
        
    return False, "No aromatic primary alcohol group found in the molecule"


# For testing purposes (can be removed/commented out when used as a module):
if __name__ == "__main__":
    test_smiles = [
        "CCCC1=C(C(O)=O)C(CO)=CO1",   # 4-(hydroxymethyl)-2-propylfuran-3-carboxylic acid (expected: True)
        "CCCCn1c(CO)ccc1C=O",         # 1-butylpyrraline (expected: True)
        "CC1=CC=C(CO)C(C)=C1",        # 2,4-dimethylbenzyl alcohol (expected: True)
        "C1=CC(=CN=C1)CO",            # 3-pyridinemethanol (expected: True)
        "CCCC",                      # simple alkane (expected: False)
    ]
    for smi in test_smiles:
        result, reason = is_aromatic_primary_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")
"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: Guaiacols
Defined as: Any phenol carrying an additional methoxy substituent at the ortho‐position.
This code scans through every aromatic carbon in the molecule (after adding hydrogens)
and checks if any bond between two aromatic carbons has one end substituted by a hydroxyl (–OH) 
and the other substituted by a methoxy (–OCH3) group.
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is defined as any phenol molecule carrying an additional methoxy group (–OCH3)
    on an ortho-position relative to the –OH.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains a guaiacol moiety, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help us count the attached hydrogens reliably
    mol = Chem.AddHs(mol)
    
    # Helper: Check if an aromatic carbon carries a hydroxyl (-OH) substituent.
    def is_hydroxy_substituted(atom):
        # We require atom to be carbon and aromatic.
        if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
            return False
        # Look at neighbors that are not aromatic carbons (i.e. substituents)
        for nbr in atom.GetNeighbors():
            # We want an oxygen that is not part of an aromatic system.
            if nbr.GetSymbol() == "O" and not nbr.GetIsAromatic():
                # For hydroxyl, the oxygen should be bonded to a hydrogen.
                # Check if at least one attached atom is H
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetSymbol() == "H":
                        # Also check that the oxygen's degree is exactly 2 (attached to the aromatic carbon and H)
                        if nbr.GetDegree() == 2:
                            return True
        return False

    # Helper: Check if an aromatic carbon carries a methoxy (-OCH3) substituent.
    def is_methoxy_substituted(atom):
        if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
            return False
        for nbr in atom.GetNeighbors():
            # Look for oxygen that is not in an aromatic ring (the substituent)
            if nbr.GetSymbol() == "O" and not nbr.GetIsAromatic():
                # The oxygen should be connected to one additional atom (the methyl carbon)
                # besides this aromatic carbon.
                substituents = [n for n in nbr.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
                if len(substituents) != 1:
                    continue
                carbon = substituents[0]
                if carbon.GetAtomicNum() != 6:
                    continue
                # Check that the carbon appears to be a methyl group (usually having 3 hydrogens attached)
                # Here we use GetTotalNumHs() which counts both explicit and implicit hydrogens.
                if carbon.GetTotalNumHs() == 3:
                    return True
        return False

    # Now, iterate over all bonds between aromatic carbons.
    # If one end is hydroxy-substituted and the other is methoxy-substituted, we have a guaiacol moiety.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Both atoms must be aromatic carbons.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and a1.GetIsAromatic() and a2.GetIsAromatic():
            # Check both possibilities.
            if is_hydroxy_substituted(a1) and is_methoxy_substituted(a2):
                return True, "Contains a guaiacol moiety (phenol with an ortho-methoxy substituent)"
            if is_methoxy_substituted(a1) and is_hydroxy_substituted(a2):
                return True, "Contains a guaiacol moiety (phenol with an ortho-methoxy substituent)"
    
    # If no such pair is found, classification fails.
    return False, "Does not contain a guaiacol moiety: missing phenol with ortho-methoxy substituent"

# Example usage (for testing)
if __name__ == "__main__":
    # A few test SMILES: guaiacol itself and one that is not.
    test_molecules = [
        ("COc1ccccc1O", "guaiacol"),               # should be True
        ("c1ccccc1O", "phenol without methoxy"),     # should be False
        ("COC(=O)\\C=C\\c1ccc(O)c(OC)c1", "trans-methylferulate")  # example from above
    ]
    for smi, name in test_molecules:
        result, reason = is_guaiacols(smi)
        print(f"SMILES: {smi}  NAME: {name}  -> {result}  REASON: {reason}")
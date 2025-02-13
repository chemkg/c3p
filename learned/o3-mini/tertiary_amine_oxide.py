"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: Tertiary Amine Oxide
A tertiary amine oxide is defined as an N-oxide where there are three organic groups bonded to the nitrogen atom.
This means that one nitrogen atom should carry a formal positive charge, be bonded to one oxygen atom carrying a negative charge 
(via a single bond) and three organic substituents. For additional robustness we require that each organic substituent’s first atom 
(which is directly bonded to the nitrogen) is either aliphatic (sp3-hybridized) or is part of an aromatic system.
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    The criteria are:
    • The molecule contains at least one nitrogen atom with a formal charge of +1.
    • That nitrogen atom is bonded to exactly four atoms.
    • One neighbor is an oxygen atom with a formal charge of -1 and is connected by a single bond.
    • The other three substituents (neighbors) must be carbon atoms.
    • Each of those carbon atoms must be either sp3-hybridized or aromatic.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains a tertiary amine oxide meeting the above criteria, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop over all atoms to find a candidate nitrogen atom.
    for atom in mol.GetAtoms():
        # Look only at nitrogen atoms.
        if atom.GetSymbol() != "N":
            continue
        # Check for the expected formal charge on the nitrogen.
        if atom.GetFormalCharge() != 1:
            continue
        # A tertiary amine oxide nitrogen should be bonded to four atoms.
        if atom.GetDegree() != 4:
            continue
        
        # Initialize flags and lists for neighbors.
        oxide_found = False
        organic_neighbors = []
        
        # Loop over neighbors of the nitrogen.
        for neigh in atom.GetNeighbors():
            # Get the bond between nitrogen and this neighbor.
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neigh.GetIdx())
            # If neighbor is oxygen with -1 charge and the bond is a single bond,
            # mark that we've found the oxidized oxygen.
            if neigh.GetSymbol() == "O" and neigh.GetFormalCharge() == -1:
                # Check that the N-O bond is a single bond.
                if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                oxide_found = True
            else:
                organic_neighbors.append(neigh)
                
        # Make sure we have exactly one oxide and the remaining three are supposed to be organic.
        if not oxide_found:
            continue
        if len(organic_neighbors) != 3:
            continue
        
        # To be considered organic in this scheme, the atom directly attached to the nitrogen must be carbon
        # AND either be sp3-hybridized (alkyl, for example) or be aromatic (aryl).
        for sub in organic_neighbors:
            if sub.GetSymbol() != "C":
                # Reject if any substituent is not carbon.
                return False, "Found an N-oxide but one substituent is not carbon-based"
            # Accept if the carbon is aromatic.
            if not sub.GetIsAromatic():
                # If not aromatic, require that it is sp3-hybridized.
                if sub.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    return False, "One substituent is not a simple organic group (must be sp3 or aromatic)"
        
        # If we have satisfied all criteria, then this molecule contains a tertiary amine oxide.
        return True, "Found tertiary amine oxide with N(+)-O(-) and three appropriate organic substituents."
    
    # No nitrogen matching the criteria was found.
    return False, "No tertiary amine oxide pattern found"

# Example usage:
# test_smiles = "C[N+](C)([O-])C"  # trimethylamine N-oxide should return True.
# result, reason = is_tertiary_amine_oxide(test_smiles)
# print(result, reason)
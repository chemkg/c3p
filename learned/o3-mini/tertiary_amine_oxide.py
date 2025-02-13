"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: Tertiary Amine Oxide
A tertiary amine oxide is defined as an N-oxide in which one nitrogen atom carries a formal +1 charge, is bonded
to one oxygen atom (with a formal –1 charge) via a single bond and to three organic substituents. For added
robustness, the first atom in each organic substituent (i.e. the atom directly attached to nitrogen) must be a carbon 
atom that is either sp3-hybridized (aliphatic) or aromatic.
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    The molecule must have at least one nitrogen atom meeting ALL of these criteria:
      • The nitrogen has a formal charge of +1.
      • Its total valence (explicit + implicit) is 4 and it has NO implicit hydrogens so that all four substituents are explicit.
      • Exactly one of its substituents is an oxygen atom (with a formal charge of –1) attached via a single bond.
      • The other three substituents are carbon atoms.
      • Each such carbon must be either aromatic or sp3-hybridized.
      • The bonds between the nitrogen and its carbon substituents must be single (or aromatic as flagged by RDKit).
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule meets the tertiary amine oxide criteria, False otherwise.
        str: A reason string explaining the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Loop over each atom looking for a candidate nitrogen.
    for atom in mol.GetAtoms():
        # Consider only nitrogen atoms.
        if atom.GetSymbol() != "N":
            continue

        # Check that nitrogen has the formal charge +1.
        if atom.GetFormalCharge() != 1:
            continue

        # Use total degree (explicit + implicit) to get the complete valence.
        total_deg = atom.GetTotalDegree()
        if total_deg != 4:
            continue
        
        # Reject if there are any implicit hydrogens; we want all substituents to be explicit.
        if atom.GetNumImplicitHs() != 0:
            continue

        # Now check the substituents.
        oxide_found = False
        organic_neighbors = []  # Will store tuples of (neighbor, bond)
        for neigh in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neigh.GetIdx())
            # Look for the oxygen substituent:
            if neigh.GetSymbol() == "O" and neigh.GetFormalCharge() == -1:
                # Require that the N-O bond is a single bond.
                if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                oxide_found = True
            else:
                organic_neighbors.append((neigh, bond))
                
        # Must have exactly one oxide and exactly three other (organic) substituents.
        if not oxide_found or len(organic_neighbors) != 3:
            continue

        # Check each organic substituent:
        valid_substituents = True
        for sub, bond in organic_neighbors:
            # The first atom must be carbon.
            if sub.GetSymbol() != "C":
                valid_substituents = False
                reason = "Found an N-oxide but one substituent is not carbon-based"
                break
            # The bond from the N to this carbon should be a single bond.
            # (Note: aromatic bonds are represented as single bonds with an aromatic flag in RDKit.)
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                if not bond.GetIsAromatic():
                    valid_substituents = False
                    reason = "One N–C bond is not a single (or aromatic) bond"
                    break
            # Check that the attached carbon is either sp3-hybridized or aromatic.
            if not sub.GetIsAromatic():
                if sub.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    valid_substituents = False
                    reason = "One substituent carbon is neither sp3-hybridized nor aromatic"
                    break
        if not valid_substituents:
            continue

        # If all criteria have been met, report a correct classification.
        return True, "Found tertiary amine oxide with N(+)-O(-) and three appropriate organic substituents."

    # If no nitrogen meets all the criteria, classification fails.
    return False, "No tertiary amine oxide pattern found"
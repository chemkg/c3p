"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: Hydroxynaphthoquinone
Definition: A naphthoquinone is a naphthalene (fused bicyclic aromatic system) in which at least two 
            C(=O) groups are attached to the ring (the oxygen atoms should lie outside of the core) 
            and a hydroxynaphthoquinone is any such system having at least one hydroxy (-OH) substituent 
            on one of the core carbons.
            
Approach:
  1. Parse the SMILES and add explicit hydrogens so that –OH groups are visible.
  2. Use a SMARTS query for a naphthalene core (c1ccc2ccccc2c1).
  3. For each atom in the found naphthalene core, count if it has a double bond to oxygen (carbonyl) 
     where the oxygen atom is not part of the core.
  4. Also, count hydroxy substituents (an oxygen attached to the core that has at least one hydrogen).
  5. Return True if there are at least two such carbonyl functions and at least one hydroxy substitution.
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines whether a given molecule is a hydroxynaphthoquinone based on its SMILES string.
    
    A hydroxynaphthoquinone is defined as a naphthoquinone (i.e. a naphthalene core that is substituted 
    by at least two carbonyl groups) with at least one hydroxy (-OH) substituent on the fused ring system.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a hydroxynaphthoquinone.
        str: Detailed explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that -OH groups become explicit.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a naphthalene core.
    # This pattern looks for a fused bicyclic aromatic system with 10 ring atoms.
    naphthalene_smarts = "c1ccc2ccccc2c1"
    naph_core_query = Chem.MolFromSmarts(naphthalene_smarts)
    if naph_core_query is None:
        return False, "Error in constructing naphthalene SMARTS"
    
    # Look for substructure matches for the naphthalene core.
    matches = mol.GetSubstructMatches(naph_core_query)
    if not matches:
        return False, "No fused aromatic ring system (naphthalene core) found"
    
    # For our purposes, we consider the first match (if more than one, each candidate could be inspected).
    core = set(matches[0])
    
    carbonyl_count = 0
    hydroxy_count = 0

    # Look at each atom in the candidate core.
    for idx in core:
        atom = mol.GetAtomWithIdx(idx)
        # We only consider aromatic carbons in the core.
        if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
            continue
        
        # Flag to indicate if this atom has a carbonyl (bonded to an oxygen via a double bond)
        has_carbonyl = False
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in core):
                    has_carbonyl = True
                    break
        if has_carbonyl:
            carbonyl_count += 1

        # Check for hydroxy substituent: an oxygen neighbor (not in the core) with at least one acidic hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in core):
                # We check if the oxygen has any attached H atoms.
                if nbr.GetTotalNumHs() > 0:
                    hydroxy_count += 1
                    break  # Count each core atom only once.
                    
    # For a naphthoquinone, we require at least two carbonyls attached to the naphthalene core.
    if carbonyl_count < 2:
        return False, "Core naphthalene system does not have at least two carbonyl groups (found {})".format(carbonyl_count)
    
    # For a hydroxynaphthoquinone, we require at least one hydroxy group on the core.
    if hydroxy_count < 1:
        return False, "Naphthoquinone core has {} carbonyl group(s) but lacks a hydroxy (-OH) substituent on the core".format(carbonyl_count)
    
    return True, ("Found naphthalene core with {} carbonyl group(s) and {} hydroxy substituent(s) on the core, "
                  "qualifying as a hydroxynaphthoquinone".format(carbonyl_count, hydroxy_count))

# Example usage:
# print(is_hydroxynaphthoquinone("Oc1ccc(O)c2C(=O)C=CC(=O)c12"))
"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: Tertiary Amine Oxide
A tertiary amine oxide is defined as an N-oxide in which one nitrogen atom (with a formal +1 charge)
is bonded to one oxygen atom (with a formal –1 charge) and three carbon atoms. Moreover, the first
atom in each organic substituent must be a carbon that is either sp3-hybridized or aromatic.
The classification uses a SMARTS substructure search to capture the desired motif.
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    
    The molecule must contain an N-oxide functional group where one nitrogen atom meets ALL of these criteria:
      • It has four explicit bonds and a formal charge of +1.
      • Exactly one of its substituents is an oxygen (with a formal charge of –1) attached via a single bond.
      • The other three substituents are carbon atoms.
      • Each such carbon substituent (i.e. the atom directly attached to the nitrogen) is either aromatic
        or sp3-hybridized.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule meets the tertiary amine oxide criteria; False otherwise.
        str: A reason string explaining the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that captures a nitrogen with 4 bonds having exactly three carbon substituents and one oxygen.
    # [N+X4] ensures the N is tetravalent and carries a formal +1 charge.
    # ([#6]) means that the substituent is any atom with atomic number 6 (carbon).
    # [O-] requires an oxygen atom with a formal charge of -1.
    pattern = Chem.MolFromSmarts("[N+X4]([#6])([#6])([#6])[O-]")
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No tertiary amine oxide pattern matched"
    
    # Check each match until one of them qualifies.
    for match in matches:
        # The match returns a tuple of atom indices corresponding to:
        # (nitrogen, carbon1, carbon2, carbon3, oxygen)
        idx_N = match[0]
        idx_c1 = match[1]
        idx_c2 = match[2]
        idx_c3 = match[3]
        idx_O = match[4]
        
        # Retrieve the nitrogen atom to be sure it has exactly 4 explicit neighbors.
        n_atom = mol.GetAtomWithIdx(idx_N)
        if n_atom.GetTotalDegree() != 4:
            continue  # defensive check
        
        # For each carbon neighbor, confirm that it is either aromatic or sp3-hybridized.
        valid = True
        for c_idx in (idx_c1, idx_c2, idx_c3):
            c_atom = mol.GetAtomWithIdx(c_idx)
            # Check that the atom is a carbon.
            if c_atom.GetAtomicNum() != 6:
                valid = False
                reason = "One substituent attached to the N is not carbon-based"
                break
            # The attached carbon must be aromatic OR sp3-hybridized.
            if not (c_atom.GetIsAromatic() or c_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3):
                valid = False
                reason = "One substituent carbon is neither aromatic nor sp3-hybridized"
                break
        if valid:
            return True, "Found tertiary amine oxide with N(+)-O(-) and three appropriate organic substituents."
    
    return False, "Pattern found but one or more substituents do not fulfill the required carbon characteristics."

# Example test invocation:
if __name__ == "__main__":
    # Testing one of the true positive examples: trimethylamine N-oxide SMILES: "C[N+](C)([O-])C"
    tp_smiles = "C[N+](C)([O-])C"
    result, reason = is_tertiary_amine_oxide(tp_smiles)
    print(f"SMILES: {tp_smiles}\nClassified: {result}\nReason: {reason}")
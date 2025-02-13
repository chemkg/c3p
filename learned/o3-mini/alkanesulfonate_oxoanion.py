"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
Definition: 'An alkanesulfonate in which the carbon at position 1 is attached to R, which can represent hydrogens, a carbon chain, or other groups.'
This improved function looks for an sp3, acyclic carbon directly attached (by a single bond) to a sulfonate group: S(=O)(=O)[O-].
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule contains an alkanesulfonate oxoanion moiety.
    The function searches for a sulfonate group (-S(=O)(=O)[O-]) that is connected via a single bond
    to an sp3 carbon. In addition, the carbon must be acyclic (i.e. not part of a ring).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a valid alkanesulfonate oxoanion is found, False otherwise.
        str: Detailed reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an sp3 (CX4) carbon directly attached to a sulfonate group S(=O)(=O)[O-].
    # The pattern "[CX4]-S(=O)(=O)[O-]" requires a single bond between the sp3 carbon and sulfur.
    pattern = Chem.MolFromSmarts("[CX4]-S(=O)(=O)[O-]")
    if pattern is None:
        return False, "Error constructing SMARTS pattern"
    
    # Find all matches to the pattern.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No substructure matching an sp3 carbon attached to S(=O)(=O)[O-] was found"
    
    valid_matches = 0
    reasons = []
    # For each match, verify that the carbon atom (the first atom in the pattern) is acyclic.
    # The match tuple gives the indices for atoms corresponding to [CX4] and S(=O)(=O)[O-] respectively.
    for match in matches:
        carbon_idx = match[0]  # first atom in our SMARTS pattern ([CX4])
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        if carbon_atom.IsInRing():
            reasons.append(f"Matched carbon atom at index {carbon_idx} is in a ring; skipping this match.")
            continue  # Exclude this match because it is not a typical open-chain (alkane) carbon.
        valid_matches += 1

    if valid_matches == 0:
        return False, "All matches found have the sp3 carbon in a ring, so no valid alkanesulfonate oxoanion was identified"
    
    return True, f"Found {valid_matches} valid alkanesulfonate oxoanion instance(s) in the molecule"

# Example test harness (uncomment the following lines to test locally):
# if __name__ == "__main__":
#     test_smiles = [
#         "C(CS([O-])(=O)=O)NC(C)=O",       # acetyltaurine(1-)
#         "OCCN(CCO)CCS([O-])(=O)=O",         # 2-[bis(2-hydroxyethyl)amino]ethanesulfonate
#         "CS([O-])(=O)=O",                  # methanesulfonate
#         "[O-]S(C[C@H](C(=O)[H])O)(=O)=O",   # L-3-sulfolactaldehyde(1-)
#         "C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCC",# N-dodecanoyltaurine(1-)
#         "[Na+].[H]C(=C([H])C([H])=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc3ccccc3c2C1(C)C)C([H])=C([H])C([H])=C1N(CCCCS([O-])(=O)=O)c2ccc3ccccc3c2C1(C)C",  # indocyanine green (contains two instances)
#         "O=C(C[S+]1CCCC1)c1ccccc1.[O-]S(=O)(=O)C(F)(F)C(=O)OC12CC3CC(C1)C(=O)C(C3)C2",  # a false positive example
#     ]
#     for smi in test_smiles:
#         result, reason = is_alkanesulfonate_oxoanion(smi)
#         print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")
"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
Definition: 'An alkanesulfonate in which the carbon at position 1 is attached to R, which can represent hydrogens, a carbon chain, or other groups.'
We define this as an sp3, acyclic carbon (with at least one hydrogen) that is directly attached via a single bond
to a sulfonate group S(=O)(=O)[O-]; additionally, the S atom should not be in a ring and should have the expected connectivity.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule contains an alkanesulfonate oxoanion moiety.
    This function looks for a sulfonate group (-S(=O)(=O)[O-]) connected to an sp3 (CX4) carbon.
    In addition, the attached sp3 carbon must be acyclic and have at least one hydrogen, and the sulfur
    atom should also not be part of a ring and have exactly four neighbors (three oxygens and one carbon).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a valid alkanesulfonate oxoanion is found, False otherwise.
        str: Detailed reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern: an sp3 (CX4) carbon connected by a single bond to sulfur in a sulfonate group.
    # Note: The pattern itself ("[CX4]-S(=O)(=O)[O-]") checks for the connection and the sulfonate.
    pattern = Chem.MolFromSmarts("[CX4]-S(=O)(=O)[O-]")
    if pattern is None:
        return False, "Error constructing SMARTS pattern"
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No substructure matching an sp3 carbon attached to S(=O)(=O)[O-] was found"
    
    valid_matches = 0
    reasons = []  # for collecting reasons why a specific match was rejected
    for match in matches:
        # match[0] = index of the sp3-carbon; match[1] = index of the sulfur atom
        carbon = mol.GetAtomWithIdx(match[0])
        sulfur = mol.GetAtomWithIdx(match[1])
        # Check that the carbon is acyclic
        if carbon.IsInRing():
            reasons.append(f"Carbon at idx {match[0]} is in a ring; skipping this match.")
            continue
        # Check that the carbon has at least one hydrogen.
        # (Using explicit hydrogens may require adding hydrogens first, so we rely on GetTotalNumHs())
        if carbon.GetTotalNumHs() < 1:
            reasons.append(f"Carbon at idx {match[0]} has no hydrogens; skipping this match.")
            continue
        # Check that the sulfur is not in a ring.
        if sulfur.IsInRing():
            reasons.append(f"Sulfur at idx {match[1]} is in a ring; skipping this match.")
            continue
        # Check that sulfur has exactly 4 neighbors.
        if sulfur.GetDegree() != 4:
            reasons.append(f"Sulfur at idx {match[1]} does not have 4 neighbors (has {sulfur.GetDegree()}); skipping.")
            continue
        valid_matches += 1

    if valid_matches == 0:
        detail = " ; ".join(reasons) if reasons else "No valid alkanesulfonate oxoanion instance found."
        return False, detail

    return True, f"Found {valid_matches} valid alkanesulfonate oxoanion instance(s) in the molecule"

# Example test harness (you can uncomment these lines to test locally):
# if __name__ == "__main__":
#     test_smiles = [
#         "C(CS([O-])(=O)=O)NC(C)=O",       # acetyltaurine(1-)
#         "OCCN(CCO)CCS([O-])(=O)=O",         # 2-[bis(2-hydroxyethyl)amino]ethanesulfonate
#         "CS([O-])(=O)=O",                  # methanesulfonate
#         "[O-]S(C[C@H](C(=O)[H])O)(=O)=O",   # L-3-sulfolactaldehyde(1-)
#         "C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCC",# N-dodecanoyltaurine(1-)
#         "[Na+].[H]C(=C([H])C([H])=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc3ccccc3c2C1(C)C)C([H])=C([H])C([H])=C1N(CCCCS([O-])(=O)=O)c2ccc3ccccc3c2C1(C)C",  # indocyanine green
#         "O=C(C[S+]1CCCC1)c1ccccc1.[O-]S(=O)(=O)C(F)(F)C(=O)OC12CC3CC(C1)C(=O)C(C3)C2",  # false positive example previously
#     ]
#     for smi in test_smiles:
#         result, reason = is_alkanesulfonate_oxoanion(smi)
#         print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")
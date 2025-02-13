"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
#!/usr/bin/env python
"""
Classifies: alkanesulfonate oxoanion
Definition: 'An alkanesulfonate in which the carbon at position 1 is attached to R, which can represent hydrogens, a carbon chain, or other groups.'
We define this as an sp3, acyclic carbon (with at least one hydrogen) that is directly attached via a single bond
to a sulfonate group S(=O)(=O)[O-]. In addition, the S atom should not be part of any ring and must have an oxygen connectivity
of two S=O (double bonds) and one S–O (single bond) where the singly-bound oxygen carries a –1 formal charge.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule contains an alkanesulfonate oxoanion group.
    
    The algorithm works as follows:
      1. Parse the SMILES and add hydrogens so that implicit H counts are correct.
      2. Look for a match to a refined SMARTS pattern capturing an sp3, acyclic carbon having at least one hydrogen
         that is single-bonded to a sulfur.
      3. For each match, check:
            - The carbon is acyclic and has at least 1 hydrogen.
            - The sulfur is acyclic and has degree 4.
            - The bond between the carbon and sulfur is single.
            - The sulfur has exactly three oxygen neighbors.
            - Of those oxygens, two must be double-bonded to sulfur and one must be single-bonded and have a –1 formal charge.
      4. Return True if at least one valid instance is found; otherwise return False with detailed reasons.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a valid alkanesulfonate oxoanion instance is found, False otherwise.
        str: Detailed reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adding hydrogens so that we can reliably count them
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern capturing an sp3, acyclic carbon that has >= 1 hydrogen,
    # bonded (by a single bond) to a sulfur (which we will further check)
    # The pattern [C;!R;X4;H1,H2,H3] means an sp3 carbon, not in a ring, with at least one hydrogen.
    # We then require a single bond (-) to [S;!R] meaning sulfur not in a ring.
    pattern = Chem.MolFromSmarts("[C;!R;X4;H1,H2,H3]-[S;!R]")
    if pattern is None:
        return False, "Error constructing SMARTS pattern"
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No substructure matching an sp3 carbon (with ≥1 H) attached to a non-ring sulfur was found"
    
    valid_matches = 0
    rejection_reasons = []
    
    # For each substructure match, perform extra connectivity checks
    for match in matches:
        carbon_idx, sulfur_idx = match[0], match[1]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        sulfur = mol.GetAtomWithIdx(sulfur_idx)
        
        # Check that the carbon is acyclic (redundant given the SMARTS, but extra-safe)
        if carbon.IsInRing():
            rejection_reasons.append(f"Carbon at idx {carbon_idx} is in a ring; skipping this match.")
            continue
        # Check that carbon has at least one hydrogen. GetTotalNumHs() automatically includes implicit H.
        if carbon.GetTotalNumHs() < 1:
            rejection_reasons.append(f"Carbon at idx {carbon_idx} has no hydrogens; skipping this match.")
            continue
        
        # Check that the bond between carbon and sulfur is a single bond.
        bond = mol.GetBondBetweenAtoms(carbon_idx, sulfur_idx)
        if bond.GetBondType() != BondType.SINGLE:
            rejection_reasons.append(f"Bond between carbon idx {carbon_idx} and sulfur idx {sulfur_idx} is not single; skipping.")
            continue
        
        # Check that sulfur is acyclic
        if sulfur.IsInRing():
            rejection_reasons.append(f"Sulfur at idx {sulfur_idx} is in a ring; skipping this match.")
            continue
        # Check that sulfur has exactly 4 neighbors (expected: one carbon and three oxygens)
        if sulfur.GetDegree() != 4:
            rejection_reasons.append(f"Sulfur at idx {sulfur_idx} does not have 4 neighbors (has {sulfur.GetDegree()}); skipping.")
            continue
        
        # Collect oxygen neighbors of sulfur
        oxygen_neighbors = [atom for atom in sulfur.GetNeighbors() if atom.GetAtomicNum() == 8]
        if len(oxygen_neighbors) != 3:
            rejection_reasons.append(f"Sulfur at idx {sulfur_idx} does not have exactly 3 oxygen neighbors (has {len(oxygen_neighbors)}); skipping.")
            continue
        
        # Now inspect the S–O bonds: we expect exactly two double bonds and one single bond.
        dbl_bond_count = 0
        sgl_bond_count = 0
        valid_single_oxygen = False  # for the oxygen with single bond, check formal charge
        for oxy in oxygen_neighbors:
            obond = mol.GetBondBetweenAtoms(sulfur_idx, oxy.GetIdx())
            if obond.GetBondType() == BondType.DOUBLE:
                dbl_bond_count += 1
            elif obond.GetBondType() == BondType.SINGLE:
                sgl_bond_count += 1
                # Check that the oxygen in the single bond carries a -1 formal charge.
                if oxy.GetFormalCharge() == -1:
                    valid_single_oxygen = True
            else:
                # Any other bond type is unexpected
                rejection_reasons.append(f"Unexpected bond type between sulfur idx {sulfur_idx} and oxygen idx {oxy.GetIdx()}; skipping match.")
                valid_single_oxygen = False
                break
        if dbl_bond_count != 2 or sgl_bond_count != 1 or not valid_single_oxygen:
            rejection_reasons.append(
                f"Sulfur at idx {sulfur_idx} does not have the expected oxygen connectivity (2 double bonds and 1 single bond with -1 charge); "
                f"found {dbl_bond_count} double bonds and {sgl_bond_count} single bonds (valid single O: {valid_single_oxygen}); skipping."
            )
            continue
        
        # If we reached here, this match fulfills all criteria.
        valid_matches += 1

    if valid_matches == 0:
        detail = " ; ".join(rejection_reasons) if rejection_reasons else "No valid alkanesulfonate oxoanion instance found."
        return False, detail

    return True, f"Found {valid_matches} valid alkanesulfonate oxoanion instance(s) in the molecule"

# Example test harness (uncomment to run locally):
# if __name__ == "__main__":
#     test_smiles = [
#         "C(CS([O-])(=O)=O)NC(C)=O",       # acetyltaurine(1-)
#         "OCCN(CCO)CCS([O-])(=O)=O",         # 2-[bis(2-hydroxyethyl)amino]ethanesulfonate
#         "CS([O-])(=O)=O",                  # methanesulfonate
#         "[O-]S(C[C@H](C(=O)[H])O)(=O)=O",   # L-3-sulfolactaldehyde(1-)
#         "C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCC",# N-dodecanoyltaurine(1-)
#         "FC(S([O-])(=O)=O)(F)F",            # triflate (should be false due to no H on C)
#     ]
#     for smi in test_smiles:
#         result, reason = is_alkanesulfonate_oxoanion(smi)
#         print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")
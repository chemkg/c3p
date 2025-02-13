"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
In these molecules the glycerol backbone has:
  - The sn-1 position: a CH2 group attached via an oxygen to an acyl (fatty acid) chain (i.e. a carbonyl-containing ester).
  - The sn-3 position: a CH2 group attached via an oxygen to a phosphoethanolamine head group.
  - The sn-2 position: a chiral center with (R)-configuration.
  
The approach:
  1. Parse the molecule and assign stereochemistry.
  2. Verify that an acyl ester (for sn-1) and a phosphate branch (for sn-3) exist.
  3. Loop over atoms to find a carbon with an assigned CIP code "R" and that has two carbon neighbors.
     For one neighbor (candidate for sn-1): check if it is bound via an oxygen (ester oxygen) to a carbonyl.
     For the other neighbor (candidate for sn-3): check if it is bound via an oxygen (ether oxygen) to a phosphorus.
  4. If such a chiral center is found, classify the molecule as a valid 1-acyl-sn-glycero-3-phosphoethanolamine.
  
If any of these do not hold, return False along with the failure reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    It verifies that the molecule contains:
      1. An acyl ester branch at the sn-1 position (the acyl chain is connected via an oxygen bonded to a carbonyl group).
      2. A phosphoethanolamine branch at the sn-3 position (the branch eventually connects to a phosphorus atom).
      3. A glycerol chiral center (sn-2) having an (R)-configuration.
      
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule fits the class and False otherwise;
                     the second element explains the classification decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Helper function: check if a given oxygen is part of an acyl ester branch.
    # We expect: CH2 -- O -- C(=O) ...
    def is_acyl_branch(branch_atom):
        # branch_atom is expected to be a carbon (CH2) that should be connected to an oxygen that in turn is connected to a carbonyl.
        for neigh in branch_atom.GetNeighbors():
            # Look for an oxygen neighbor (ester oxygen)
            if neigh.GetAtomicNum() == 8:
                # Check neighbors of this oxygen excluding branch_atom
                for n2 in neigh.GetNeighbors():
                    if n2.GetIdx() == branch_atom.GetIdx():
                        continue
                    # n2 should be a carbon that is part of a carbonyl, i.e. bonded double to an oxygen.
                    if n2.GetAtomicNum() == 6:
                        bond = mol.GetBondBetweenAtoms(neigh.GetIdx(), n2.GetIdx())
                        if bond is None:
                            continue
                        if bond.GetBondTypeAsDouble() == 2:
                            return True
        return False

    # Helper function: check if a given branch (a carbon) leads to the phosphate branch.
    # We expect: CH2 -- ? -- O -- ... eventually attached to a phosphorus atom.
    def is_phosphate_branch(branch_atom):
        for neigh in branch_atom.GetNeighbors():
            # Sometimes the branch_atom might be attached to an oxygen directly
            if neigh.GetAtomicNum() == 8:
                # Check if any neighbor of this oxygen (other than branch_atom) is a phosphorus.
                for n2 in neigh.GetNeighbors():
                    if n2.GetIdx() == branch_atom.GetIdx():
                        continue
                    if n2.GetAtomicNum() == 15:  # phosphorus
                        return True
        return False

    # Now search for a candidate glycerol chiral center.
    valid_chiral_found = False
    for atom in mol.GetAtoms():
        # We are interested in carbon atoms that have an assigned stereochemistry.
        if atom.GetAtomicNum() != 6 or not atom.HasProp('_CIPCode'):
            continue
        cip = atom.GetProp('_CIPCode')
        # Only consider if the CIP code is "R"
        if cip != "R":
            continue

        # A glycerol's central (sn-2) carbon should be connected to three neighbors.
        # In our expected molecules, two of these neighbors are CH2 groups (sn-1 and sn-3).
        neighbor_carbons = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                neighbor_carbons.append(nbr)
        if len(neighbor_carbons) < 2:
            continue  # not enough carbon neighbors to be glycerol backbone

        # Among the carbon neighbors, one should have an acyl ester branch and one should be connected to a phosphate.
        found_acyl = False
        found_phosphate = False
        for nbr in neighbor_carbons:
            # Ensure that nbr is a CH2; glycerol sn-1 and sn-3 are typically CH2 groups.
            # (Using degree >=2 as a rough filter)
            if nbr.GetDegree() < 2:
                continue
            if not found_acyl and is_acyl_branch(nbr):
                found_acyl = True
            if not found_phosphate and is_phosphate_branch(nbr):
                found_phosphate = True
        
        if found_acyl and found_phosphate:
            valid_chiral_found = True
            break

    if not valid_chiral_found:
        return False, "No glycerol chiral center with (R)-configuration having both acyl and phosphate branches found"

    return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine with (R)-configuration"

# Example usage
if __name__ == "__main__":
    # Test with one example: 1-stearoyl-sn-glycero-3-phosphoethanolamine
    test_smiles = "CCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"
    result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(test_smiles)
    print(result, reason)
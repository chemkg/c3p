"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
In these molecules the glycerol backbone has:
  - sn-1: a CH2 group linked via an oxygen to an acyl (fatty acid) chain (i.e. through an ester-connected carbonyl)
  - sn-3: a CH2 group linked via an oxygen to a phosphoethanolamine headgroup
  - sn-2: a chiral center that is required to have (R)-configuration.
  
The approach:
  1. Parse the molecule and assign stereochemistry.
  2. Identify candidate chiral carbons (with explicit chirality) that are connected to two carbon neighbors (representing sn-1 and sn-3).
  3. For one neighbor, verify that an acyl branch exists (by finding an oxygen which in turn is attached to a carbonyl group).
  4. For the other neighbor, verify that a phosphate branch exists (by finding an oxygen that eventually leads to a phosphorus atom).
  5. If such a candidate is found, check its assigned CIP code; report failure if it is not “R”.
  
If the steps do not lead to a valid candidate we return (False, <reason>).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    It verifies that the molecule contains:
      1. A glycerol backbone in which the sn-1 position (one CH2 group) carries an acyl ester branch 
         (i.e. an oxygen that connects to a carbonyl group).
      2. The sn-3 position (the other CH2 group) connects via oxygen to a phosphate (P).
      3. A chiral center at sn-2 which should have (R)-configuration based on CIP rules.
      
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule fits the class, False otherwise;
                     the second element provides the reason for the classification.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Clean and assign stereochemistry (force re-assignment to get CIP codes)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Helper function: check if a given carbon (candidate branch of glycerol) leads to an acyl ester branch.
    # We expect that branch_atom (a CH2) is attached to an oxygen that is itself bound (double-bonded) to a carbonyl.
    def is_acyl_branch(branch_atom):
        for neigh in branch_atom.GetNeighbors():
            # Look for an oxygen (ester oxygen)
            if neigh.GetAtomicNum() == 8:
                for n2 in neigh.GetNeighbors():
                    if n2.GetIdx() == branch_atom.GetIdx():
                        continue
                    # n2 should be a carbon that is connected via a double bond to an oxygen.
                    if n2.GetAtomicNum() == 6:
                        bond = mol.GetBondBetweenAtoms(neigh.GetIdx(), n2.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() == 2:
                            return True
        return False

    # Helper function: check if a given branch (a CH2 group) leads to a phosphate branch.
    # We expect it to have an oxygen neighbor that is in turn connected to a phosphorus atom.
    def is_phosphate_branch(branch_atom):
        for neigh in branch_atom.GetNeighbors():
            if neigh.GetAtomicNum() == 8:
                for n2 in neigh.GetNeighbors():
                    if n2.GetIdx() == branch_atom.GetIdx():
                        continue
                    if n2.GetAtomicNum() == 15:  # phosphorus
                        return True
        return False

    # Loop over atoms to find a candidate glycerol chiral center.
    # For the glycerol backbone in these molecules, the chiral center (sn-2) is a carbon
    # that is bonded (as heavy atoms) to one oxygen (hydroxyl) and two carbons (the sn-1 and sn-3 CH2 groups).
    candidate_found = False
    candidate_reason = ""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # We require an explicit chirality flag on the central carbon.
        if atom.GetChiralTag() == Chem.CHI_UNSPECIFIED:
            continue
        # Try to get the CIP code if available.
        cip = ""
        if atom.HasProp('_CIPCode'):
            cip = atom.GetProp('_CIPCode')
        else:
            cip = "unknown"
        
        # Among heavy atom neighbors, collect carbon neighbors (sn-1 and sn-3 should be carbons).
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # We expect exactly two carbon heavy neighbors (sn-1 and sn-3) for a glycerol backbone center.
        if len(carbon_neighbors) != 2:
            continue

        # For each of the two carbon neighbors, test for acyl ester and phosphate branches.
        found_acyl = False
        found_phosphate = False
        for nbr in carbon_neighbors:
            # We use degree (number of connections) as a rough filter.
            if nbr.GetDegree() < 2:
                continue
            # Check if this branch leads to an acyl chain.
            if not found_acyl and is_acyl_branch(nbr):
                found_acyl = True
                continue  # move to next neighbor
            # Check if this branch leads to a phosphate branch.
            if not found_phosphate and is_phosphate_branch(nbr):
                found_phosphate = True

        if found_acyl and found_phosphate:
            candidate_found = True
            # Now check the CIP configuration.
            if cip != "R":
                return False, f"Glycerol chiral center found but configuration is {cip} instead of R"
            else:
                return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine with (R)-configuration"
        else:
            # Record details for debugging if needed.
            candidate_reason = "Candidate chiral center failed: " + \
                               f"acyl branch found: {found_acyl}, phosphate branch found: {found_phosphate}"
    
    if not candidate_found:
        return False, "No glycerol chiral center with both acyl and phosphate branches found. " + candidate_reason

    return False, "Molecule did not match criteria for 1-acyl-sn-glycero-3-phosphoethanolamine"

# Example usage (testing one sample):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"  # 1-stearoyl-sn-glycero-3-phosphoethanolamine
    result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(test_smiles)
    print(result, reason)
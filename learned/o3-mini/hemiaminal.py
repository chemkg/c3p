"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemiaminal
Definition: “Any organic amino compound that has an amino group and a hydroxy group attached 
to the same carbon atom. Hemiaminals are intermediates in the formation of imines by addition of 
an amine to an aldehyde or ketone; those derived from primary amines are particularly unstable.”
 
This version improves upon the previous one by enforcing that the nitrogen atom attached to the candidate 
carbon must have at least one hydrogen (using its total hydrogen count) in addition to not being double-bonded 
to oxygen. This helps to remove many false positives.
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has a tetrahedral (mostly sp³) carbon that bears both a hydroxyl group (-OH)
    and an amino group (-NH2, -NHR, or -NR2) attached directly to that same carbon.
    
    The strategy is:
      1. Parse the molecule and add explicit hydrogens.
      2. Iterate over sp³ carbons with exactly 4 bonds.
      3. For each carbon, check that at least one neighbor is a genuine –OH:
           - oxygen must have exactly 2 neighbors (the candidate carbon and one hydrogen)
      4. Also, check that at least one neighbor is an amino nitrogen:
           - the bond to the candidate carbon is single,
           - the nitrogen is not double‐bonded to an oxygen (indicative of an amide),
           - and the nitrogen has at least one hydrogen (from explicit+implicit count).
    
    If any candidate carbon meets both criteria, the molecule is declared to possess a hemiaminal group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a hemiaminal is detected, False otherwise.
        str: Reason/message with matching atom indices or explanation.
    """
    # Parse SMILES and check validity.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string; could not parse molecule"
        
    # Add explicit hydrogens so that –OH groups and NH counts are more reliably detected.
    mol = Chem.AddHs(mol)
    
    hemiaminal_matches = []  # will collect tuples (carbon_idx, oxygen_idx, nitrogen_idx)

    # Helper: Check if an oxygen is a bona fide hydroxyl (-OH)
    def is_hydroxyl(oxygen, parent_idx):
        # In a free hydroxyl, oxygen should have exactly two neighbors: the candidate carbon and one hydrogen.
        nbrs = oxygen.GetNeighbors()
        if len(nbrs) != 2:
            return False
        valid_neighbor = False
        has_h = False
        for nbr in nbrs:
            if nbr.GetAtomicNum() == 1:
                has_h = True
            else:
                # The heavy atom neighbor must be the candidate carbon.
                if nbr.GetIdx() != parent_idx:
                    return False
                else:
                    valid_neighbor = True
        return has_h and valid_neighbor

    # Helper: Check if a nitrogen is a free amino group (non-amide, and with at least one hydrogen)
    def is_valid_amino(nitrogen):
        # First, ensure that nitrogen has at least one hydrogen (implicit or explicit)
        if nitrogen.GetTotalNumHs() < 1:
            return False
        # Now check that there is no double bond from nitrogen to an oxygen (which may indicate an amide)
        for nbr in nitrogen.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), nbr.GetIdx())
            if bond is None:
                continue
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                return False
        return True

    # Examine every carbon in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # only consider carbon atoms
        # Require sp³ hybridization.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # With all hydrogens explicit, a tetrahedral carbon should have exactly 4 neighbors.
        if atom.GetDegree() != 4:
            continue

        current_idx = atom.GetIdx()
        hydroxyl_matches = []  # oxygen neighbors valid as –OH
        amino_matches = []     # nitrogen neighbors valid as amino groups
        
        for nbr in atom.GetNeighbors():
            # Consider only single bonds.
            bond = mol.GetBondBetweenAtoms(current_idx, nbr.GetIdx())
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            atomic_num = nbr.GetAtomicNum()
            if atomic_num == 8:
                # Check oxygen neighbor for valid hydroxyl group.
                if is_hydroxyl(nbr, current_idx):
                    hydroxyl_matches.append(nbr.GetIdx())
            elif atomic_num == 7:
                # Check nitrogen neighbor for free amino character.
                if is_valid_amino(nbr):
                    amino_matches.append(nbr.GetIdx())
        
        # If at least one hydroxyl and one amino neighbor are found, count this carbon as a hemiaminal center.
        if hydroxyl_matches and amino_matches:
            # Record all combinations.
            for o_idx in hydroxyl_matches:
                for n_idx in amino_matches:
                    hemiaminal_matches.append((current_idx, o_idx, n_idx))
    
    if hemiaminal_matches:
        return True, f"Found hemiaminal substructure in atoms with indices: {hemiaminal_matches}"
    else:
        return False, "No carbon with both hydroxyl and amino substituents found; not a hemiaminal"

# Example usage (testing a couple of examples):
if __name__ == "__main__":
    test_smiles_list = [
        "OC(N)CC",  # 2-Aminopropanol: expected to be hemiaminal.
        "C1[C@@]2(N3CC[C@@]42[C@]5(N(C=6C4=CC=CC6)C(C[C@]7([C@@]5([C@@]1(C(=CCO7)C3)[H])[H])[H])=O)[H])O"  # pseudostrychnine: a known hemiaminal example.
    ]
    for s in test_smiles_list:
        result, reason = is_hemiaminal(s)
        print(result, reason)
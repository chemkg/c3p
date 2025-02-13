"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoid – Members of the class of tetrahydrochromenochromene.
Definition: Rotenoids consist of a cis‐fused tetrahydrochromenochromene skeleton and its substituted derivatives.
Our improved approach uses refined criteria based on ring analysis. We look for a fused bicyclic system defined by:
  1. Two rings that share at least two atoms (i.e. a shared bond).
  2. One ring that is highly aromatic (chromene-like) and one that is partially saturated.
     We quantify aromaticity by the fraction of atoms in the ring that are declared aromatic.
  3. At least one of the rings must contain an oxygen atom.
This refined heuristic should reduce false classifications.
"""

from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is likely to be a rotenoid based on its SMILES string.
    It uses an improved heuristic:
      1. Retrieve all rings (as lists of atom indices).
      2. Check pairs of rings that are fused (share at least 2 atoms = a whole bond).
      3. Compute the aromatic fraction for each ring.
         A ring is considered mostly aromatic if (aromatic atoms/size) >= 0.75.
         A ring is considered partially saturated if its aromatic fraction is lower.
      4. If one ring meets the aromatic criteria and the other does not,
         and if at least one of the two rings contains an oxygen atom,
         then classify the molecule as a rotenoid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule likely belongs to the rotenoid class, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information: each ring is a tuple of atom indices.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop over all pairs of rings looking for a fused bicyclic system.
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            # For rings to be fused as rings (sharing at least one bond), they should share at least 2 atoms.
            shared_atoms = set(rings[i]).intersection(rings[j])
            if len(shared_atoms) < 2:
                continue  # Not fused strongly enough
            
            # Gather atom objects for each ring.
            ring1_atoms = [mol.GetAtomWithIdx(idx) for idx in rings[i]]
            ring2_atoms = [mol.GetAtomWithIdx(idx) for idx in rings[j]]
            
            # Calculate aromatic fraction for each ring.
            arom1 = sum(1 for atom in ring1_atoms if atom.GetIsAromatic())
            arom2 = sum(1 for atom in ring2_atoms if atom.GetIsAromatic())
            frac1 = arom1 / len(rings[i])
            frac2 = arom2 / len(rings[j])
            
            # Define thresholds:
            # A "chromene-like" ring should be mostly aromatic (at least 75% aromatic atoms).
            # The partner (tetrahydro-like) ring should not reach this threshold.
            aromatic_threshold = 0.75
            
            cond1 = (frac1 >= aromatic_threshold and frac2 < aromatic_threshold)
            cond2 = (frac2 >= aromatic_threshold and frac1 < aromatic_threshold)
            
            if not (cond1 or cond2):
                continue  # This pair does not show the desired difference in aromaticity
            
            # Check that at least one of the rings contains an oxygen atom.
            has_oxygen_ring1 = any(atom.GetAtomicNum() == 8 for atom in ring1_atoms)
            has_oxygen_ring2 = any(atom.GetAtomicNum() == 8 for atom in ring2_atoms)
            if not (has_oxygen_ring1 or has_oxygen_ring2):
                continue  # Neither ring contains oxygen
            
            # If the pair meets all criteria, we classify as rotenoid.
            return True, ("Molecule contains a fused bicyclic system "
                          "with one predominantly aromatic ring and one partially saturated ring, "
                          "with at least one ring bearing oxygen, consistent with a rotenoid skeleton")
    
    # If no pair was found, return false.
    return False, "Molecule does not contain the fused bicyclic system characteristic of rotenoids"

# Example usage:
if __name__ == "__main__":
    # Test examples (only a few given here; more can be added for further validation)
    test_examples = {
        "Millettosin": "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C6=C(OC(C=C6)(C)C)C=C5",
        "LSM-4404": "CC(C)C1CC2=C(O1)C=CC3=C2OC4=C(C3=O)C5=CC(=C(C=C5OC4)OC)OC",
        "tephrosin": "[H][C@@]12COc3cc(OC)c(OC)cc3[C@]1(O)C(=O)c1ccc3OC(C)(C)C=Cc3c1O2",
        "Stemonone (FN)": "O1C2=C(C=3C(OC2=O)=CC(OC)=C(OC)C3)C(=O)C=4C1=CC(OC)=CC4O"
    }
    for name, smi in test_examples.items():
        result, reason = is_rotenoid(smi)
        print(f"{name}: {result} -- {reason}")
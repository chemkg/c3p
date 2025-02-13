"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin pigments – a subclass of flavonoids
Anthoxanthins are typically flavones or flavonols that contain a 2-phenylchromen-4-one core.
They have a fused system (rings A and C with a keto group at position 4) with an attached benzene ring (ring B),
often decorated with hydroxyl and/or methoxy groups. They are water‐soluble pigments (i.e. low lipophilicity)
with colors ranging from colorless/white to yellow.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen

def is_anthoxanthin(smiles: str):
    """
    Determines whether the molecule given as a SMILES string is likely to be an anthoxanthin.
    
    We use the following criteria:
      1. Valid SMILES string.
      2. The molecule must contain a flavonoid core defined as a 2-phenylchromen-4-one scaffold.
         In our approach the flavonoid core is identified by a SMARTS pattern for the fused (rings A and C)
         part that contains the keto group.
      3. There must be an aromatic benzene ring (ring B) directly attached to the core.
      4. There must be at least one free hydroxyl group.
      5. The Crippen logP must be low (<= 3.0) in keeping with water‐solubility.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as an anthoxanthin,
                     and False otherwise. The second element is the reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more specific SMARTS pattern for the flavonoid core
    # This pattern is meant to capture the fused benzopyranone system (rings A and C, with a C4 keto group)
    # Note: substituents are not specified so extra groups are tolerated.
    flavonoid_core = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(c2)")
    core_matches = mol.GetSubstructMatches(flavonoid_core)
    if not core_matches:
        return False, "No flavonoid core (2-phenylchromen-4-one) found"
    # Take the first match (assuming one core per molecule)
    core_atoms = set(core_matches[0])
    
    # Now, look for an attached aromatic benzene ring.
    # Get all rings: we consider rings that are aromatic, contain exactly 6 atoms, and (if possible) are pure carbon rings.
    ring_info = mol.GetRingInfo()
    aromatic_benzene_rings = []
    for ring in ring_info.AtomRings():
        # Check if it is aromatic and has 6 atoms
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Check if all atoms in the ring are carbon (atomic number 6)
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                aromatic_benzene_rings.append(set(ring))
    
    # An attached benzene ring (ring B) should not be completely part of the core.
    # We require that at least one of its atoms is bonded to an atom in the flavonoid core.
    attached_benzene_count = 0
    for ring in aromatic_benzene_rings:
        if ring.issubset(core_atoms):
            continue  # skip rings fully inside the core match
        # For each atom in the ring, check bonds to see if any neighbor is in the core.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms and nbr.GetIdx() not in ring:
                    attached_benzene_count += 1
                    break
            else:
                continue
            break  # found an attached ring, move to next ring

    if attached_benzene_count != 1:
        return False, f"Expected exactly one benzene ring attached to the flavonoid core, found {attached_benzene_count}"
    
    # Check for free hydroxyl groups (–OH)
    hydroxyl = Chem.MolFromSmarts("[OX2H]")
    n_hydrox = len(mol.GetSubstructMatches(hydroxyl))
    if n_hydrox < 1:
        return False, "Lacks free hydroxyl group(s) typical of anthoxanthins"
    
    # Calculate Crippen logP to assess water-solubility.
    # Anthoxanthins are water-soluble pigments; we require a logP <= 3.0.
    logP = Crippen.MolLogP(mol)
    if logP > 3.0:
        return False, f"LogP too high ({logP:.2f}); molecule appears too lipophilic for a typical anthoxanthin"
    
    return True, "Molecule contains a specific flavonoid core with one attached benzene ring, hydroxylation, and low logP, consistent with an anthoxanthin"

# Example usage (uncomment to test):
# print(is_anthoxanthin("COC1=CC=C(C=C1O)C1=CC(=O)C2=C(O)C=C(O)C=C2O1"))  # diosmetin
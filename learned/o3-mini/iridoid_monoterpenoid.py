"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
#!/usr/bin/env python
"""
Classifies: Iridoid Monoterpenoid
Definition: Iridoid monoterpenoids are biosynthesized from isoprene and usually have a 
cyclopentane ring fused to a six-membered oxygen heterocycle. (Secoiridoids arise when 
one bond in the cyclopentane ring is cleaved.) 

This program uses ring analysis rather than pure SMARTS matching to identify a fused 
bicyclic core consisting of one 5-membered and one 6-membered ring that share at least 
two atoms. In the 6-membered ring, at least one oxygen atom must be present.
In addition, the total number of carbons in the fused system is expected to lie in a 
narrow range (7–9 carbons). This aims to reduce both false positives (when other parts 
of a larger molecule trigger a match) and false negatives (where the core is decorated).
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines whether the molecule (SMILES string) is an iridoid monoterpenoid.
    The strategy:
      1. Parse the SMILES into a molecule.
      2. Get ring information from the molecule.
      3. Look for fused rings where one is 5-membered and one is 6-membered with at least two atoms in common.
      4. Require that the 6-membered ring contains at least one oxygen.
      5. Take the union of the atoms in both rings and count the carbon atoms.
      6. Return True only if the number of carbons in the bicyclic core is between 7 and 9.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as an
                     iridoid monoterpenoid and False otherwise; the second element is a message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Search for a pair of fused rings: one 5-membered and one 6-membered with overlap.
    for i in range(len(rings)):
        ring1 = rings[i]
        for j in range(i+1, len(rings)):
            ring2 = rings[j]
            # Check for one ring of size 5 and one of size 6.
            if (len(ring1) == 5 and len(ring2) == 6) or (len(ring1) == 6 and len(ring2) == 5):
                # They must share at least two atoms to be considered fused.
                shared_atoms = set(ring1) & set(ring2)
                if len(shared_atoms) < 2:
                    continue

                # Identify the 6-membered ring for oxygen-check.
                ring6 = ring1 if len(ring1) == 6 else ring2
                # Check that the six-membered ring contains at least one oxygen.
                oxygen_count = sum(1 for idx in ring6 if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                if oxygen_count < 1:
                    continue

                # Now consider the fused system as the union of both rings.
                fused_atoms = set(ring1) | set(ring2)
                # Count the number of carbon atoms in the fused system.
                carbon_count = sum(1 for idx in fused_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if 7 <= carbon_count <= 9:
                    return True, (f"Found fused bicyclic core (5-membered and 6-membered with oxygen) "
                                  f"with {carbon_count} carbons")
                else:
                    # If a fused core is found but the carbon count is off, then reject.
                    return False, (f"Fused bicyclic core found but carbon count is outside expected range "
                                   f"({carbon_count} carbons found)")
    
    return False, "No fused bicyclic iridoid monoterpenoid core was detected"


# (Optional) Test cases using example SMILES from the prompt:
if __name__ == "__main__":
    test_cases = [
        ("Hamigeran L", "BrC1=C(O)C(=C([C@@H]2[C@@](CC[C@@H]2C(C)C)(CC(=O)O)C)C=C1C)C(=O)O"),
        ("Lippioside I", "OC1(C2C(CC1O)C(=COC2OC3OC(C(O)C(O)C3O)COC(=O)/C=C/C4=CC=C(O)C=C4)C(O)=O"),
        ("methyl 7alpha-acetoxydeacetylbotryoloate", "O=C(OC)[C@@H]1[C@]2(O)[C@]([C@H](OC(=O)C)C([C@@H]2[C@@H](O)C[C@H]1C)(C)C)(CO)C"),
        ("Valechlorin", "ClCC1(O)C2C(=CC1OC(=O)CC(C)C)C(=COC2OC(=O)CC(C)C)COC(=O)C"),
        ("Shanzhiside methyl ester", "O[C@@]1([C@@]2([C@]([C@H](O)C1)(C(=CO[C@H]2O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)C(OC)=O)[H])[H])C")
    ]
    
    for name, smi in test_cases:
        result, reason = is_iridoid_monoterpenoid(smi)
        print(f"Name: {name}")
        print("SMILES:", smi)
        print("Result:", result, "| Reason:", reason)
        print("-" * 80)
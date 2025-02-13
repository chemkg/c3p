"""
Classifies: CHEBI:28802 flavonols
"""
#!/usr/bin/env python
"""
Classifies: Flavonols 
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic ring is replaced by a hydroxy group.
In our approach we first check for a flavone (2-phenylchromen-4-one) core, and then check that in the pyran ring
there is at least one carbon (other than the carbonyl carbon) that carries an external hydroxyl (–OH) group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is defined as a hydroxyflavone in which the ring hydrogen at position 3
    of the heterocyclic (pyran) ring is replaced by a hydroxy group.
    
    Our strategy:
      1. Parse the SMILES.
      2. Check that the molecule contains a flavone core. We use a SMARTS
         pattern corresponding to the basic 2-phenylchromen-4-one scaffold.
         (This is a heuristic check and may not capture all cases.)
      3. Examine the rings in the molecule. Among rings of size 6, we select one that:
           - Contains one oxygen (the heterocyclic O)
           - Contains one carbonyl carbon (a carbon that is double-bonded to O, but that O is exocyclic)
         This ring is taken as the candidate pyran ring.
      4. Finally, in this candidate ring we look for at least one carbon (other than the carbonyl) 
         that has an –OH substituent attached (and the hydroxyl is not part of the ring).
         This extra –OH is assumed to be in the 3‑position.
         
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a flavonol, False otherwise.
        str: Reason detailing the result.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a flavone (2-phenylchromen-4-one) core.
    # This pattern matches an aromatic benzene ring attached to a heterocycle containing an oxygen and a carbonyl.
    flavone_smarts = Chem.MolFromSmarts("c1ccccc1-c2oc3ccccc3c(=O)c2")
    if not mol.HasSubstructMatch(flavone_smarts):
        return False, "Molecule does not contain the basic flavone (2-phenylchromen-4-one) core"

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    # Flag to indicate if we find a pyran ring (6-membered ring with one oxygen and one carbonyl)
    found_pyran_with_OH = False

    # Iterate over all rings
    for ring in ring_info:
        if len(ring) != 6:
            continue  # we are looking for a 6-membered ring

        num_ring_oxygen = 0
        carbonyl_indices = []  # indices of carbons in the ring that are carbonyl carbons

        # Check each atom in the ring for oxygen and carbonyl properties
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "O":
                num_ring_oxygen += 1
            # Identify a carbon that is double-bonded to an oxygen (carbonyl)
            if atom.GetAtomicNum() == 6:
                for nb in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                    # Check if bond is a double bond and neighbor is an oxygen
                    if nb.GetAtomicNum() == 8 and bond.GetBondType().name == "DOUBLE":
                        # To be safe, require that the oxygen is not part of the ring
                        if nb.GetIdx() not in ring:
                            carbonyl_indices.append(atom.GetIdx())
                            break

        # Heuristic: the candidate pyran ring should have exactly one ring oxygen and one carbonyl carbon.
        if num_ring_oxygen == 1 and len(carbonyl_indices) == 1:
            # Now search for a hydroxyl (OH) substituent attached to any carbon in this ring that is not the carbonyl.
            for idx in ring:
                if idx in carbonyl_indices:
                    continue  # skip the carbonyl carbon (position 4)
                atom = mol.GetAtomWithIdx(idx)
                # We are interested in aromatic carbons (sp2) that could be at position 3.
                if atom.GetAtomicNum() == 6:
                    # Look at neighbors of this carbon (exclude atoms that are part of the ring)
                    for nb in atom.GetNeighbors():
                        # Check if neighbor is an oxygen (and not in the ring) with at least one attached hydrogen (i.e. –OH)
                        if nb.GetAtomicNum() == 8 and nb.GetIdx() not in ring and nb.GetTotalNumHs() > 0:
                            # Found a candidate: an external hydroxyl substituent on a ring carbon.
                            found_pyran_with_OH = True
                            break
                if found_pyran_with_OH:
                    break
        if found_pyran_with_OH:
            break

    if not found_pyran_with_OH:
        return False, "Flavone core detected but lacks the required 3-hydroxy substitution on the heterocyclic ring"

    return True, "Molecule contains a flavonol skeleton with a 3-hydroxy substitution on the heterocyclic ring"


# Example usage (you can remove or comment these out in production):
if __name__ == "__main__":
    # Test a couple of example SMILES strings.
    test_smiles = [
        # kaempferol SMILES from the examples
        "Oc1ccc(cc1)-c1oc2cc(O)cc(O)c2c(=O)c1O",
        # flavone (which is not flavonol because the 3-position is not hydroxy)
        "c1ccccc1-c2oc3ccccc3c(=O)c2"
    ]
    for smi in test_smiles:
        result, reason = is_flavonols(smi)
        print(f"SMILES: {smi}")
        print(f"Flavonol? {result} ({reason})\n")
"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone family – esterified prostanoids from marine corals.
Definition: A class of esterified prostanoids obtained from marine corals.
Heuristic criteria:
  • Must contain a cyclopentenone ring (i.e. a five‐membered ring containing at least one C=O group 
    where the oxygen is exocyclic, and at least one carbon–carbon double bond in the ring).
  • Must contain at least one ester group (–O–C(=O)–).
Note: The patterns below are a heuristic approximation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone (esterified prostanoid from marine corals)
    based on its SMILES string.

    The function uses heuristic criteria:
      1. The SMILES can be parsed.
      2. The molecule contains a cyclopentenone ring, defined as any five-membered ring 
         that includes at least one ketone group (a carbon double-bonded to oxygen, where the O is not part of the ring)
         and at least one carbon–carbon double bond within the ring.
      3. The molecule contains at least one ester functional group.
      4. Optionally checks that the molecular weight is in a plausible range (>300 Da).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a clavulone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, search for a cyclopentenone-type ring.
    # We iterate over all five-membered rings and check for two features:
    #   (a) A ketone group: a carbon in the ring that has a double bond to an oxygen atom 
    #       that is NOT part of the ring (i.e. the O is exocyclic).
    #   (b) At least one carbon–carbon double bond in the ring (C=C).
    ring_info = mol.GetRingInfo().AtomRings()
    cyclopentenone_found = False
    for ring in ring_info:
        if len(ring) != 5:
            continue  # only consider 5-membered rings
        ketone = False
        cc_double = False

        # Check for ketone: for every atom in the ring, if carbon then look for an exocyclic C=O.
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                # Check bond is a double bond
                if bond.GetBondTypeAsDouble() == 2.0:
                    # Identify the neighbor atom which is not the current atom.
                    nbr = bond.GetOtherAtom(atom)
                    # Check if neighbor is oxygen and NOT in the ring.
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                        ketone = True
                        break
            if ketone:
                break

        # Check for at least one carbon-carbon double bond within the ring.
        # We iterate over every pair of distinct atoms in the ring.
        for i in range(len(ring)):
            for j in range(i+1, len(ring)):
                bond = mol.GetBondBetweenAtoms(ring[i], ring[j])
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    # Ensure that both atoms are carbons (i.e., we are not counting the ketone bond)
                    a1 = mol.GetAtomWithIdx(ring[i])
                    a2 = mol.GetAtomWithIdx(ring[j])
                    if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                        cc_double = True
                        break
            if cc_double:
                break

        if ketone and cc_double:
            cyclopentenone_found = True
            break

    if not cyclopentenone_found:
         return False, "No cyclopentenone ring found (required prostanoid core missing)"

    # Next, check for at least one ester group via SMARTS.
    ester_smarts = "[OX2][CX3](=O)[#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if not mol.HasSubstructMatch(ester_pattern):
         return False, "No ester group found (molecule is not esterified)"

    # Optional check: ensure the molecular weight is plausible for clavulones (>300 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
         return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical clavulone"

    return True, "Molecule contains a cyclopentenone ring and at least one ester group consistent with clavulones"

# Example usage:
if __name__ == "__main__":
    # Test example using one of the provided SMILES.
    test_smiles = "BrC1=C[C@@](O)([C@@]([C@@H](OC(=O)C)/C=C/CCCC(OC)=O)(C1=O)[H])C/C=C\\CCCCC"
    result, reason = is_clavulone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
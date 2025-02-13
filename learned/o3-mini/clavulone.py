"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone family – esterified prostanoids from marine corals.
Definition: A class of esterified prostanoids obtained from marine corals.
Heuristic criteria (improved):
  • Molecule must be parseable from SMILES.
  • It must contain a non-aromatic cyclopentenone ring. We search for a five-membered ring that:
       – Has at least one carbon which is double-bonded to oxygen (ketone functionality). 
         (We allow the oxygen to be exocyclic to the ring.)
       – Contains at least one carbon–carbon double bond between ring atoms.
  • It must contain at least one ester group (–O–C(=O)–).
  • To reduce false positives we also require that the total number of rings is not too high (≤ 3).
Note: This heuristic will not be perfect, and some true clavulones may be missed or some false positives remain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone (esterified prostanoid from marine corals)
    based on its SMILES string and heuristic structural features.
    
    The criteria used herein are:
      1. SMILES can be parsed into a molecule.
      2. The molecule contains a non-aromatic cyclopentenone ring. This ring is defined as:
           - A five-membered ring.
           - At least one carbon in the ring that has a double bond to oxygen (ketone),
             where the oxygen is not itself in the ring.
           - At least one carbon–carbon double bond within the ring.
           - None of the ring atoms are aromatic.
      3. The molecule contains at least one ester group (–O–C(=O)–).
      4. The total number of rings in the molecule is not too high (≤ 3).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a clavulone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get rings from the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    cyclopentenone_found = False
    for ring in ring_info:
        if len(ring) != 5:
            continue  # only consider five-membered rings
        
        # Check that none of the ring atoms are aromatic
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        ketone = False
        cc_double = False
        
        # Check for ketone functionality: at least one carbon in the ring must have a double bond to oxygen,
        # where that oxygen is not in the ring.
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue  # we expect a carbon
            for bond in atom.GetBonds():
                # Check for double bond
                if bond.GetBondTypeAsDouble() == 2.0:
                    nbr = bond.GetOtherAtom(atom)
                    # Accept if neighbor is oxygen and is not in the ring.
                    if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in ring):
                        ketone = True
                        break
            if ketone:
                break
        
        # Check for at least one carbon–carbon double bond within the ring.
        # Iterate over pairs of ring atom indices.
        n = len(ring)
        for i in range(n):
            for j in range(i+1, n):
                bond = mol.GetBondBetweenAtoms(ring[i], ring[j])
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    a1 = mol.GetAtomWithIdx(ring[i])
                    a2 = mol.GetAtomWithIdx(ring[j])
                    # Ensure that both atoms are carbons
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
    
    # Check that the molecule contains at least one ester group.
    # This SMARTS pattern matches an oxygen single-bonded to a carbon that is double-bonded to oxygen,
    # with the carbon attached to at least one other carbon.
    ester_smarts = "[OX2][CX3](=O)[#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if not mol.HasSubstructMatch(ester_pattern):
         return False, "No ester group found (molecule is not esterified)"
    
    # To reduce false positives, check that the total number of rings is not too high.
    total_rings = mol.GetRingInfo().NumRings()
    if total_rings > 3:
         return False, f"Too many rings ({total_rings}) for a typical clavulone"

    return True, "Molecule contains a non-aromatic cyclopentenone ring and at least one ester group consistent with clavulones"

# Example usage:
if __name__ == "__main__":
    # Test a known example (7-Acetoxy-7,8-dihydrobromovulone II)
    test_smiles = "BrC1=C[C@@](O)([C@@]([C@@H](OC(=O)C)/C=C/CCCC(OC)=O)(C1=O)[H])C/C=C\\CCCCC"
    result, reason = is_clavulone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
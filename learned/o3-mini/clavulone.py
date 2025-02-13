"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone family – esterified prostanoids from marine corals.
Definition: A class of esterified prostanoids obtained from marine corals.
Heuristic criteria (improved):
  • Molecule must be parseable from SMILES.
  • The overall molecular weight should be relatively high (≥300 Da)
    and the molecule should have at least 15 carbon atoms.
  • It must contain a non‐aromatic cyclopentenone ring.
      – If the five‐membered ring has an in‐ring carbonyl conjugated to a C=C bond
        (or if the carbonyl is exocyclic), we regard that as the “enone” prostanoid core.
      – For epoxidized derivatives (where the double bond is replaced by an epoxide),
        we relax the requirement on the C=C bond if a ketone functionality is present.
  • It must contain at least one ester group (–O–C(=O)–).
  • To reduce false positives we also require that the total number of rings is not too high (≤ 3).
Note: This heuristic won’t perfectly capture every clavulone, but should help improve the F1 score.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone (esterified prostanoid from marine corals)
    based on its SMILES string and heuristic structural features.
    
    Criteria:
      1. SMILES must be parsed.
      2. The molecule should be large enough: MW >=300 Da and at least 15 carbon atoms.
      3. The molecule should contain a non‐aromatic cyclopentenone ring. We try two SMARTS:
             a. One that finds a five‐membered ring containing an in‐ring carbonyl and a C=C bond.
             b. One that finds a five‐membered ring with a substituent carbonyl (exocyclic).
         If neither matches, we (if an epoxide group is present anywhere)
         iterate over five‐membered non‐aromatic rings and look for at least one ketone bond.
      4. The molecule must have at least one ester group.
      5. The total number of rings is limited (≤ 3).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a clavulone, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check global size: molecular weight and carbon count.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a clavulone"
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 15:
        return False, f"Too few carbons ({n_carbons}) for a typical clavulone"
    
    # Prepare SMARTS patterns for cyclopentenone:
    # Pattern 1: five‐membered non‐aromatic ring with an in‐ring ketone and a C=C bond.
    cp_pattern1 = Chem.MolFromSmarts("[C;r5;!a]1[C;r5;!a]=[C;r5;!a][C;r5;!a](=O)[C;r5;!a]1")
    # Pattern 2: five‐membered non‐aromatic ring where the carbonyl is exocyclic
    cp_pattern2 = Chem.MolFromSmarts("[C;r5;!a]1[C;r5;!a]=[C;r5;!a]([C;r5;!a]1)C(=O)")
    
    cyclopentenone_found = False
    if mol.HasSubstructMatch(cp_pattern1) or mol.HasSubstructMatch(cp_pattern2):
        cyclopentenone_found = True
    else:
        # Sometimes the enone is altered by epoxidation.
        # If an epoxide is present anywhere (SMARTS for epoxide ring: three‐membered ring with oxygen),
        epoxide_pattern = Chem.MolFromSmarts("[OX2r3]")
        has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
        # Then attempt a manual search: iterate over every five‐membered non‐aromatic ring,
        # and check if at least one carbon in the ring has a double bond to oxygen (ketone)
        ring_info = mol.GetRingInfo().AtomRings()
        for ring in ring_info:
            if len(ring) != 5:
                continue
            if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            ketone = False
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Only carbons can form the needed ketone
                if atom.GetAtomicNum() != 6:
                    continue
                # Look for a double bond to oxygen where the oxygen is not part of the ring.
                for bond in atom.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2.0:
                        nbr = bond.GetOtherAtom(atom)
                        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                            ketone = True
                            break
                if ketone:
                    break
            if ketone and has_epoxide:
                cyclopentenone_found = True
                break

    if not cyclopentenone_found:
         return False, "No cyclopentenone ring found (required prostanoid core missing)"
    
    # Check for at least one ester group (pattern: oxygen single-bonded to carbon that is double-bonded to oxygen).
    ester_smarts = "[OX2][CX3](=O)[#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if not mol.HasSubstructMatch(ester_pattern):
         return False, "No ester group found (molecule is not esterified)"
    
    # Limit the total number of rings to reduce false positives.
    total_rings = mol.GetRingInfo().NumRings()
    if total_rings > 3:
         return False, f"Too many rings ({total_rings}) for a typical clavulone"
    
    return True, "Molecule contains a valid cyclopentenone core and ester group consistent with clavulones"

# Example usage (you may add more tests):
if __name__ == "__main__":
    # This is one example (7-Acetoxy-7,8-dihydrobromovulone II)
    test_smiles = "BrC1=C[C@@](O)([C@@]([C@@H](OC(=O)C)/C=C/CCCC(OC)=O)(C1=O)[H])C/C=C\\CCCCC"
    result, reason = is_clavulone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
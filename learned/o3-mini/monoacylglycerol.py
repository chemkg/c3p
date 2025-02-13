"""
Classifies: CHEBI:17408 monoacylglycerol
"""
#!/usr/bin/env python
"""
Classifies: monoacylglycerol

Definition:
  A monoacylglycerol is defined as a glyceride in which one of the substituents on the glycerol
  backbone is an acyl group (linked via an ester –O–C(=O)– bond) while the two remaining substituents
  are free (typically -OH groups).

Heuristic Overview:
  1. Parse the SMILES string.
  2. Look for a three‐carbon (sp³) chain that is a candidate glycerol backbone.
  3. For each candidate chain, collect all oxygen atoms bonded to its carbons that are not part of the chain.
  4. Classify these oxygens as ester oxygen (if the oxygen is bonded to a carbon that is also double bonded to another oxygen)
     or as free hydroxyl (if it carries at least one hydrogen).
  5. If exactly three oxygen substituents are found – one ester and two free – we classify the molecule as a monoacylglycerol.
  
Note:
  This is a heuristic approach. In cases where the substituents or ordering differ significantly, the heuristic
  may fail. In that case the function returns False with an explanation.
"""

from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol should contain a glycerol backbone (a linear chain of three sp³ carbons) 
    that bears three oxygen substituents. One of these oxygens should form an ester bond (attached to a 
    carbon that carries a carbonyl group) and the other two should appear as free hydroxyl groups.
    
    Args:
        smiles (str): A SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a monoacylglycerol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper function: returns True if the oxygen atom is part of an ester bond.
    def is_ester_oxygen(oxy):
        """
        Check if the oxygen 'oxy' is directly bonded to a carbon that is also double-bonded to another oxygen.
        This is taken to indicate an ester linkage (–O–C(=O)–).
        """
        # For each neighbor of the oxygen:
        for nbr in oxy.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                # Examine bonds from that carbon.
                for bond in nbr.GetBonds():
                    # Skip the bond to oxy.
                    if bond.GetOtherAtom(nbr) == oxy:
                        continue
                    # If bond type is double and the other atom is oxygen, count as carbonyl.
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(nbr).GetAtomicNum() == 8:
                        return True
        return False

    # Helper function: returns True if the oxygen appears to be a free hydroxyl.
    def is_free_hydroxyl(oxy):
        """
        Check if the oxygen atom 'oxy' appears to be a free hydroxyl group.
        We check if it carries at least one hydrogen (explicit or implicit).
        """
        if oxy.GetTotalNumHs() > 0:
            return True
        return False

    # --- Step 1: Identify candidate glycerol backbones.
    # Instead of using a strict SMARTS match (which in our previous attempt required
    # each backbone carbon to have exactly one oxygen neighbor), we can search for any three-carbon
    # chain (sp³) and then check the substituents.
    candidate_indices = []
    # Loop over all atoms; only consider carbons with sp3 hybridization.
    atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    n_atoms = len(atoms)
    # Create a mapping from atom index to atom for convenience.
    atom_idx_set = set(a.GetIdx() for a in atoms)
    
    # We now look for connected triplets a-b-c where a, b, c are distinct and connected in sequence.
    for a in atoms:
        for b in a.GetNeighbors():
            if b.GetIdx() not in atom_idx_set:
                continue
            # Ensure a and b are connected by a single (non-aromatic) bond.
            bond_ab = mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
            if bond_ab is None or bond_ab.GetBondType() != Chem.BondType.SINGLE:
                continue
            for c in b.GetNeighbors():
                if c.GetIdx() not in atom_idx_set or c.GetIdx() == a.GetIdx():
                    continue
                bond_bc = mol.GetBondBetweenAtoms(b.GetIdx(), c.GetIdx())
                if bond_bc is None or bond_bc.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # We require that a and c are not directly bonded (to ensure a linear chain of three atoms).
                if mol.GetBondBetweenAtoms(a.GetIdx(), c.GetIdx()):
                    continue
                # We have a candidate chain: record the ordered triple (a, b, c)
                candidate = (a.GetIdx(), b.GetIdx(), c.GetIdx())
                candidate_indices.append(candidate)
    
    if not candidate_indices:
        return False, "No three-carbon sp³ backbone found"

    # --- Step 2: For each glycerol candidate, check the oxygen substituents.
    for cand in candidate_indices:
        # Get the three backbone atoms.
        backbone_atoms = [mol.GetAtomWithIdx(idx) for idx in cand]
        # Collect unique oxygen atoms (neighbors not in the backbone).
        oxy_set = {}
        for atom in backbone_atoms:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in cand:
                    oxy_set[nbr.GetIdx()] = nbr
        # A proper glycerol backbone (or its modified derivative) should have exactly 3 oxygen substituents.
        if len(oxy_set) != 3:
            continue

        ester_count = 0
        free_oh_count = 0
        for oxy in oxy_set.values():
            if is_ester_oxygen(oxy):
                ester_count += 1
            elif is_free_hydroxyl(oxy):
                free_oh_count += 1
            else:
                # If the oxygen is neither clearly ester nor free hydroxyl,
                # we do not count it.
                pass

        if ester_count == 1 and free_oh_count == 2:
            return True, "Glycerol backbone found with one acyl (ester) substituent and two free hydroxyl groups"

    return False, "No glycerol backbone with exactly one acyl (ester) substitution and two free hydroxyl groups found"


# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "C(CCCCCCCCCCCC(OCC(CO)O)=O)CCCCCCC",  # 1-icosanoylglycerol
        "CCCCCC\\C=C/CCCCCCCC(=O)OC(CO)CO",     # 2-linoleoylglycerol
        "CCCCCCCCC(=O)OC[C@@H](O)CO",           # 1-decanoyl-sn-glycerol
        "O(C[C@@H](O)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",  # MG(22:5(...))
        "O=C1C2=C(O)C(=C(C)C=C2C(=O)C=3C1=C(O)C=C(O)C3)C(=O)OC[C@H](O)CO",  # A false positive candidate
    ]
    for s in test_smiles:
        res, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\n Result: {res}\n Reason: {reason}\n")
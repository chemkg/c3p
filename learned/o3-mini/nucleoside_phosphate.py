"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
#!/usr/bin/env python
"""
Classifies: nucleoside phosphate
A nucleoside phosphate is defined as:
  "A nucleobase-containing molecular entity that is a nucleoside in which one or more of the 
   sugar hydroxy groups has been converted into a mono- or poly-phosphate. The term includes 
   both nucleotides and non-nucleotide nucleoside phosphates."

This improved approach checks for three key features:
  1. A candidate sugar ring (we look for 5-membered rings with four carbons and one oxygen).
  2. A phosphate group (P atom) attached via an oxygen to the sugar ring (allowing for exocyclic attachment).
  3. A nucleobase substructure (using a set of SMARTS patterns for common purine and pyrimidine cores).
Due to the structural diversity of these molecules, this classification remains heuristic.
"""

from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate has a sugar (typically a furanose ring), a nucleobase and
    one or more phosphate group(s) attached to the sugar via an oxygen bridge.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a nucleoside phosphate, False otherwise.
        str: Explanation of the reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Identify candidate sugar rings.
    # We look for 5-membered rings (furanose) that contain exactly one oxygen
    # and four carbons. (This is a heuristic for ribose/deoxyribose)
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # list of sets of atom indices corresponding to candidate sugar rings
    for ring in ring_info:
        if len(ring) == 5:
            o_count = 0
            c_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    o_count += 1
                if atom.GetAtomicNum() == 6:
                    c_count += 1
            # Typical furanose ring: 4 carbons and 1 oxygen
            if o_count == 1 and c_count == 4:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No candidate sugar ring (5-membered ring with 4 carbons and 1 oxygen) detected"

    # 2. Check for phosphate attachment.
    # For each phosphorus atom, see if one of its neighboring oxygen atoms 
    # is connected (as a neighbor) to an atom that is part of any candidate sugar ring.
    phosphate_attached = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus atom
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                    # Check its neighbors (excluding the phosphorus itself)
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == atom.GetIdx():
                            continue
                        # If this neighbor belongs to any sugar ring, we consider the phosphate as attached
                        for sugar in sugar_rings:
                            if nbr2.GetIdx() in sugar:
                                phosphate_attached = True
                                break
                        if phosphate_attached:
                            break
                    if phosphate_attached:
                        break
            if phosphate_attached:
                break
    if not phosphate_attached:
        return False, "No phosphate group found attached to the sugar moiety"

    # 3. Check for the presence of a nucleobase.
    # We use a set of SMARTS patterns for purine and pyrimidine rings.
    nucleobase_smarts_list = [
        "n1cnc2ncnc12",       # common purine motif (adenine/guanine core)
        "c1nc2c(n1)cnc2",      # alternative purine pattern
        "n1cnc(=O)n1",        # cytosine-like pattern
        "O=c1[nH]c(=O)[nH]1",  # uracil-like pattern
        "CC1=CN(C(=O)NC1=O)"   # thymine (includes methyl) 
    ]
    nucleobase_found = False
    for smarts in nucleobase_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase substructure found (expected purine or pyrimidine ring)"

    # All three features were found.
    return True, "Contains a nucleobase, a candidate sugar (furanose) ring, and a phosphate group attached via an oxygen bridge"

# Example usage:
if __name__ == "__main__":
    # Test with a known nucleoside phosphate: UDP
    udp_smiles = "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O"
    result, reason = is_nucleoside_phosphate(udp_smiles)
    print("UDP:", result, reason)
    
    # Additional test: GTP
    gtp_smiles = "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O"
    result, reason = is_nucleoside_phosphate(gtp_smiles)
    print("GTP:", result, reason)
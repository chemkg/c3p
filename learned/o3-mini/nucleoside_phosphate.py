"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
A nucleoside phosphate is defined as:
  "A nucleobase-containing molecular entity that is a nucleoside in which one or more of the 
   sugar hydroxy groups has been converted into a mono- or poly-phosphate. The term includes 
   both nucleotides and non-nucleotide nucleoside phosphates."
   
This program checks for three features:
  1. Presence of a nucleobase (using several SMARTS patterns for common nucleobases)
  2. Presence of a sugar moiety (a five-membered ring containing exactly one oxygen)
  3. Attachment of at least one phosphate group to the sugar via an oxygen bond

Due to the structural diversity of these molecules, the classification is heuristic.
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate has a nucleobase linked to a sugar ring (typically a furanose) 
    with one or more sugar hydroxy groups converted into mono- or poly-phosphates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Find a candidate sugar ring: We look for a 5-membered ring that contains exactly one oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []
    for ring in ring_info:
        if len(ring) == 5:
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            # Most furanose sugars have exactly one ring oxygen (and four carbons)
            if o_count == 1:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No candidate sugar ring (5-membered ring with one oxygen) detected"
    
    # 2. Check for phosphate attachment to the sugar.
    # We look for at least one phosphorus atom (atomic number 15) bonded (via an intervening oxygen)
    # to any atom of a sugar ring. It is common for the phosphate to be attached via an O atom.
    phospho_attached = False
    for sugar in sugar_rings:
        # For each atom in the sugar ring, examine neighbors that are not in the ring.
        for idx in sugar:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in sugar:
                    continue
                # We expect an oxygen substituent (hydroxyl group) on the sugar
                if nbr.GetAtomicNum() == 8:
                    # Now check if this oxygen is bonded to a phosphorus.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 15:
                            phospho_attached = True
                            break
                    if phospho_attached:
                        break
            if phospho_attached:
                break
        if phospho_attached:
            break
    if not phospho_attached:
        return False, "No phosphate group found attached to the sugar moiety"
    
    # 3. Check for the presence of a nucleobase.
    # Here we use several SMARTS patterns common to typical nucleobases.
    nucleobase_smarts_list = [
        # Purine pattern (matches e.g., adenine and guanine cores)
        "n1cnc2ncnc12",
        # Cytosine
        "n1cnc(=O)n1",
        # Uracil (and can match thymine minus a methyl)
        "O=c1[nH]c(=O)[nH]1",
        # Thymine (including the methyl group)
        "CC1=CN(C(=O)NC1=O)"
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
        return False, "No nucleobase substructure found (expected purine/pyrimidine ring)"
    
    return True, "Contains a nucleobase, a candidate sugar ring, and a phosphate group attached to the sugar"

# Example usage:
if __name__ == "__main__":
    # Test with a known nucleoside phosphate: UDP
    udp_smiles = "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O"
    result, reason = is_nucleoside_phosphate(udp_smiles)
    print("UDP:", result, reason)
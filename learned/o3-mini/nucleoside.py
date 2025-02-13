"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside (an N‐glycosyl compound containing both a nucleobase and a ribose/deoxyribose sugar).

Definition:
“A nucleoside is an N-glycosyl compound that has both a nucleobase – normally adenine, guanine, xanthine, thymine, cytosine or uracil – 
and either a ribose or deoxyribose as functional parents.”
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside must contain a sugar moiety (ribose or deoxyribose) and a nucleobase 
    that are connected via an N-glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------------------------------
    # Step 1. Identify candidate sugar rings.
    # For ribose/deoxyribose the sugar is a five-membered (furanose) ring typically with 4 carbons and 1 oxygen,
    # and the atoms are saturated (sp3 hybridized).
    sugar_rings = []
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        if len(ring) == 5:
            oxy_count = 0
            carbon_count = 0
            sp3_flag = True
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxy_count += 1
                if atom.GetAtomicNum() == 6:
                    carbon_count += 1
                # Check if atom is sp3 hybridized (typical for sugars)
                if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    sp3_flag = False
            if oxy_count == 1 and carbon_count == 4 and sp3_flag:
                sugar_rings.append(ring)
    
    if not sugar_rings:
        return False, "No sugar moiety (furanose ring with 4 carbons and 1 oxygen) detected"
    
    # -------------------------------
    # Step 2. Identify the nucleobase.
    # We use a set of SMARTS patterns that are meant to capture common purine and pyrimidine rings.
    nucleobase_patterns = [
        "n1cnc2ncnc2n1",      # Generic purine pattern (e.g., adenine and guanine cores)
        "c1[nH]c(=O)nc1",      # Pattern that may capture uracil or cytosine cores
        "c1nc(=O)n(C)c1=O"     # Example pattern for thymine-like structures
    ]
    nucleobase_atoms_set = set()
    found_nucleobase = False
    for smarts in nucleobase_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                found_nucleobase = True
                for match in matches:
                    for idx in match:
                        nucleobase_atoms_set.add(idx)
    if not found_nucleobase:
        return False, "No nucleobase detected via common nucleobase patterns"
    
    # -------------------------------
    # Step 3. Check for N-glycosidic connectivity.
    # We search for a bond between an atom of the sugar ring (typically a carbon, e.g. the anomeric carbon)
    # and a nitrogen atom that is part of the nucleobase.
    for ring in sugar_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only carbon atoms from the sugar ring as the likely anomeric center.
            if atom.GetAtomicNum() != 6:
                continue
            # Check each neighbor of the sugar carbon that is not part of the sugar itself.
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in ring:
                    # Check if neighbor is a nitrogen (atomic number 7) and is part of the nucleobase match.
                    if nbr.GetAtomicNum() == 7 and nbr_idx in nucleobase_atoms_set:
                        return True, "Molecule contains an N-glycosidic bond connecting a sugar moiety to a nucleobase"
    
    return False, "Found sugar and nucleobase motifs, but no N-glycosidic bond connecting them"
"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition:
“A nucleoside is an N-glycosyl compound that has both a nucleobase – normally adenine, guanine, xanthine,
thymine, cytosine or uracil – and either a ribose or deoxyribose as functional parents.”
Our strategy:
  1. Parse the SMILES string.
  2. Exclude molecules containing phosphorus (indicative of a nucleotide rather than a nucleoside).
  3. Identify candidate sugar rings: five-membered (furanose) rings with exactly 4 carbons, 1 oxygen,
     and all atoms sp3.
  4. Identify candidate nucleobase rings: rings (of size ≥5) that are aromatic with at least two nitrogen atoms.
  5. Check for connectivity between a sugar carbon (typically the anomeric center) and a nitrogen atom of the nucleobase.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside must contain a ribose/deoxyribose sugar moiety and a nucleobase connected via an N-glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules with phosphorus (usually nucleotides with phosphate groups)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Phosphate group detected; molecule is likely a nucleotide, not a nucleoside"

    # -------------------------------
    # Step 1. Identify candidate sugar rings.
    # For a ribose or deoxyribose, expect a five-membered (furanose) ring with 4 carbons and 1 oxygen, all sp3.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []
    for ring in ring_info:
        if len(ring) == 5:
            oxy_count = 0
            carbon_count = 0
            sp3_flag = True
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxy_count += 1
                elif atom.GetAtomicNum() == 6:
                    carbon_count += 1
                else:
                    # Unexpected atom type for a sugar ring (e.g. heteroatom not typical in ribose)
                    sp3_flag = False
                if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    sp3_flag = False
            if oxy_count == 1 and carbon_count == 4 and sp3_flag:
                # record as a set for easier membership testing
                sugar_rings.append(set(ring))
                
    if not sugar_rings:
        return False, "No sugar moiety (furanose ring with 4 carbons and 1 oxygen) detected"
    
    # -------------------------------
    # Step 2. Identify candidate nucleobase rings.
    # We look for rings (size >= 5) that are aromatic and have at least two nitrogen atoms.
    nucleobase_atoms = set()
    for ring in ring_info:
        if len(ring) >= 5:
            nitrogen_count = 0
            aromatic_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    nitrogen_count += 1
                if atom.GetIsAromatic():
                    aromatic_count += 1
            # Heuristic: nucleobases are usually aromatic heterocycles.
            if nitrogen_count >= 2 and aromatic_count >= (len(ring) // 2):
                nucleobase_atoms.update(ring)
                
    if not nucleobase_atoms:
        return False, "No nucleobase detected via aromatic heterocycle patterns"
    
    # -------------------------------
    # Step 3. Check for N-glycosidic connectivity.
    # Look for a bond between a sugar ring carbon (likely the anomeric carbon) and a nitrogen atom in the nucleobase.
    for sugar in sugar_rings:
        for idx in sugar:
            atom = mol.GetAtomWithIdx(idx)
            # Focus on carbon atoms in the sugar as potential anomeric centers.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # The connecting bond should leave the sugar (nbr not in sugar ring) and land in the nucleobase.
                if nbr_idx not in sugar and nbr_idx in nucleobase_atoms:
                    # Ensure that the neighbor is nitrogen and the bond is a single bond (typical of glycosidic bonds)
                    if nbr.GetAtomicNum() == 7:
                        bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                        if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            return True, "Molecule contains an N-glycosidic bond connecting a sugar moiety to a nucleobase"
                            
    return False, "Found sugar and nucleobase motifs but no connecting N-glycosidic bond detected"
"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue – a molecule that can substitute for a normal nucleobase in nucleic acids.
Heuristic criteria used (refined):
  1. Molecule is valid and its molecular weight is in the 80–350 Da range.
  2. Contains at least 2 nitrogen atoms overall.
  3. Contains at least one aromatic heterocycle that is nucleobase‐like. Specifically,
     a) Either the molecule has a single six‐membered aromatic ring with at least 2 nitrogen atoms
        (pyrimidine‐like), OR
     b) It has a fused ring system (two rings sharing at least 2 atoms) where one ring is six‐membered
        and the other is five‐membered and the combined (unique) nitrogen count is at least 3 (purine‐like).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
    
    Heuristic criteria (refined):
      - Must parse to a valid molecule.
      - Molecular weight between 80 and 350 Da.
      - Contains at least 2 nitrogen atoms.
      - Contains at least one nucleobase-like aromatic heterocycle:
            Either a six-membered aromatic ring with at least 2 nitrogens (pyrimidine-like)
            OR a fused ring system (sharing >=2 atoms) composed of a six-membered and a five-membered ring
            with combined nitrogen count >= 3 (purine-like).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule fits the refined criteria, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is not in the typical range for nucleobase analogues (80-350 Da)"
    
    # Count total nitrogen atoms in the molecule
    total_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if total_nitrogens < 2:
        return False, f"Found only {total_nitrogens} nitrogen atom(s); nucleobase analogues typically have 2 or more"
    
    # Retrieve ring information and gather aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = []
    for ring in ring_info.AtomRings():
        # check if all atoms in the ring are aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            ring_data = {"indices": set(ring), "size": len(ring), "n_count": n_in_ring}
            aromatic_rings.append(ring_data)
    
    if not aromatic_rings:
        return False, "Molecule has no aromatic rings; nucleobase analogues are heterocyclic"
    
    candidate_found = False
    candidate_reason = ""
    
    # First, check for a single six-membered aromatic ring with at least two nitrogens (pyrimidine-like)
    for ring in aromatic_rings:
        if ring["size"] == 6 and ring["n_count"] >= 2:
            candidate_found = True
            candidate_reason = "Contains a six-membered aromatic heterocycle with at least two nitrogen atoms (pyrimidine-like)"
            break
    
    # Next, if no single ring candidate found, look for a fused ring system (purine-like)
    if not candidate_found:
        n_arom = len(aromatic_rings)
        for i in range(n_arom):
            for j in range(i+1, n_arom):
                # Check if rings share at least 2 atoms (fused)
                intersection = aromatic_rings[i]["indices"].intersection(aromatic_rings[j]["indices"])
                if len(intersection) >= 2:
                    # Look for a pair where one ring is six-membered and the other is five-membered
                    size_i = aromatic_rings[i]["size"]
                    size_j = aromatic_rings[j]["size"]
                    if ( (size_i == 6 and size_j == 5) or (size_i == 5 and size_j == 6) ):
                        # Count unique nitrogens in the union of the two rings
                        union_indices = aromatic_rings[i]["indices"].union(aromatic_rings[j]["indices"])
                        union_nitrogens = sum(1 for idx in union_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                        if union_nitrogens >= 3:
                            candidate_found = True
                            candidate_reason = "Contains a fused aromatic ring system (6- and 5-membered) with sufficient nitrogen atoms (purine-like)"
                            break
            if candidate_found:
                break
    
    if not candidate_found:
        return False, "Molecule does not contain a nucleobase-like aromatic heterocycle (neither pyrimidine-like nor purine-like)"
    
    return True, "Molecule fits refined heuristic criteria for a nucleobase analogue"
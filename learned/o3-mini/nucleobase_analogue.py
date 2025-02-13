"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue – a molecule that can substitute for a normal nucleobase in nucleic acids.
Improved heuristic criteria:
  1. Must parse as a valid molecule.
  2. Molecular weight between 80 and 350 Da.
  3. Contains at least 2 nitrogen atoms.
  4. Contains a single nucleobase-like aromatic heterocycle defined as either:
     a) A six-membered aromatic ring with at least 2 nitrogen atoms (pyrimidine‐like),
         AND no additional nitrogen-containing aromatic rings in the molecule,
         OR
     b) A fused ring system (two rings that share at least 2 atoms) where one ring is six-membered 
         and the other is five-membered and the union of the two contains at least 3 nitrogen atoms (purine‐like).
  5. The candidate aromatic system must represent a significant portion (>=50%) of the heavy atoms 
     in the molecule (to reduce misclassification of molecules with extra aromatic substituents, sugars, etc.).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    The heuristic criteria are:
      - Valid molecule and molecular weight between 80 and 350 Da.
      - Contains at least 2 nitrogen atoms.
      - Contains one nucleobase-like aromatic heterocycle, defined as:
           * Either a six-membered aromatic ring that contains ≥2 nitrogen atoms (pyrimidine-like)
           * OR a fused system of a six-membered ring and a five-membered ring sharing at least 2 atoms,
             with a union of at least 3 nitrogen atoms (purine-like).
      - The candidate heterocycle must represent at least 50% of the molecule’s heavy atoms.
      - The molecule should not contain additional unrelated nitrogen-containing aromatic systems.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets these criteria, False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (using exact weight)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is not in the typical range (80-350 Da) for nucleobase analogues"
    
    # Count total nitrogen atoms
    total_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if total_nitrogens < 2:
        return False, f"Found only {total_nitrogens} nitrogen atom(s); expected at least 2"
    
    # Get aromatic ring information from the molecule.
    ring_info = mol.GetRingInfo()
    aromatic_rings = []
    # We identify rings that are aromatic (all atoms aromatic) and record:
    # - set of atom indices in the ring
    # - ring size
    # - number of nitrogen atoms in the ring
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            aromatic_rings.append({
                "indices": set(ring),
                "size": len(ring),
                "n_count": n_in_ring
            })
    
    if not aromatic_rings:
        return False, "Molecule has no aromatic rings; nucleobase analogues are heterocyclic"
    
    candidate_type = None  # will be "pyrimidine" or "purine"
    candidate_atoms = set()  # indices for the candidate nucleobase heterocycle
    
    # First, search for a pyrimidine-like candidate: one single six-membered ring with at least 2 nitrogens.
    pyrimidine_candidates = [ring for ring in aromatic_rings if ring["size"] == 6 and ring["n_count"] >= 2]
    if len(pyrimidine_candidates) == 1:
        candidate_type = "pyrimidine-like"
        candidate_atoms = pyrimidine_candidates[0]["indices"]
    else:
        # If not found, try to find a purine-like candidate: a pair of fused rings (one 6-membered and one 5-membered)
        # sharing at least 2 atoms and with the union having at least 3 nitrogen atoms.
        n_arom = len(aromatic_rings)
        found_purine = False
        for i in range(n_arom):
            for j in range(i+1, n_arom):
                ring1 = aromatic_rings[i]
                ring2 = aromatic_rings[j]
                # Check that one ring is 6-membered and the other 5-membered.
                if not (((ring1["size"] == 6 and ring2["size"] == 5) or (ring1["size"] == 5 and ring2["size"] == 6))):
                    continue
                # For a fused system, the rings must share at least 2 atoms.
                if len(ring1["indices"].intersection(ring2["indices"])) < 2:
                    continue
                union_indices = ring1["indices"].union(ring2["indices"])
                union_nitrogens = sum(1 for idx in union_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                if union_nitrogens >= 3:
                    candidate_type = "purine-like"
                    candidate_atoms = union_indices
                    found_purine = True
                    break
            if found_purine:
                break
    
    if candidate_type is None:
        return False, "Molecule does not contain a nucleobase-like aromatic heterocycle (neither pyrimidine-like nor purine-like)"
    
    # Now, ensure that the molecule does not contain extra independent heterocyclic aromatic rings.
    # We consider only rings that contain at least one nitrogen.
    nucleo_rings = [ring for ring in aromatic_rings if ring["n_count"] > 0]
    # For a pyrimidine-like candidate, we expect only one such ring.
    if candidate_type == "pyrimidine-like" and len(nucleo_rings) > 1:
        return False, "Molecule contains extra nitrogen-containing aromatic rings beyond the pyrimidine-like core"
    # For a purine-like candidate, we expect exactly two rings in the fused system.
    if candidate_type == "purine-like" and len(nucleo_rings) > 2:
        return False, "Molecule contains extra nitrogen-containing aromatic rings beyond the purine-like fused core"
    
    # Next, check what fraction of the heavy atoms are in the candidate nucleobase heterocycle.
    total_heavy = mol.GetNumHeavyAtoms()
    frac = len(candidate_atoms) / total_heavy
    if frac < 0.5:
        return False, f"Candidate heterocycle accounts for only {frac:.2f} fraction of heavy atoms (<0.5), suggesting additional substituents"
    
    return True, f"Molecule fits refined heuristic criteria as a nucleobase analogue ({candidate_type})."
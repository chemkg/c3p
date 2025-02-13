"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue – a molecule that can substitute for a normal nucleobase in nucleic acids.
The heuristic criteria here are:
  1. Must parse as a valid molecule; molecular weight between 80 and 350 Da.
  2. Contains at least 2 nitrogen atoms.
  3. Contains one nucleobase-like aromatic heterocycle that is either:
       a) A single six-membered aromatic ring with ≥2 nitrogens (pyrimidine-like)
          (and no extra fused aromatic rings in that core),
       OR
       b) A fused system of exactly two rings – one six-membered and one five-membered that share at least 2 atoms,
          whose union contains at least 3 nitrogen atoms (purine-like).
  4. The candidate core must represent a significant fraction (>=35% of heavy atoms; relaxed from 50%) of the molecule.
  
Extra aromatic rings that are not fused to this core should not be present.
This code clusters all aromatic rings that share atoms and then tests each fused system candidate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the nucleobase analogue heuristic criteria, False otherwise.
        str: Explanation of the decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Molecular weight check (exact weight in Daltons)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is not in the 80-350 Da range typical for nucleobase analogues"
    
    # Must have at least 2 nitrogen atoms
    total_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if total_nitrogens < 2:
        return False, f"Found only {total_nitrogens} nitrogen atom(s); need at least 2"
    
    # Get ring info: we focus on aromatic rings (all atoms aromatic)
    ring_info = mol.GetRingInfo().AtomRings()
    rings = []
    for ring in ring_info:
        # Only consider rings where every atom is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            rings.append({
                "indices": set(ring),
                "size": len(ring),
                "n_count": n_in_ring
            })
    
    if not rings:
        return False, "Molecule has no aromatic rings; nucleobase analogues must be heterocyclic molecules"
    
    # Cluster rings that are fused (i.e. share at least one atom) into connected components.
    # We use a simple union–find like procedure.
    n = len(rings)
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx
    
    # For every pair of rings, if they share at least one atom, union their indices.
    for i in range(n):
        for j in range(i+1, n):
            if rings[i]["indices"].intersection(rings[j]["indices"]):
                union(i, j)
    
    # Group rings by their root parent.
    comp_dict = {}
    for i in range(n):
        r = find(i)
        comp_dict.setdefault(r, []).append(i)
        
    # List to store candidate fused systems (one candidate per connected component).
    candidates = []
    for comp in comp_dict.values():
        # For each component, combine the rings.
        comp_rings = [rings[i] for i in comp]
        # Compute union of atom indices for the fused system candidate:
        union_atoms = set()
        for ring in comp_rings:
            union_atoms = union_atoms.union(ring["indices"])
        # Sum nitrogen counts within the union (note: an atom shared between rings is counted once).
        union_nitrogens = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        # Also record the sizes of the individual rings.
        ring_sizes = [ring["size"] for ring in comp_rings]
        
        candidates.append({
            "ring_count": len(comp_rings),
            "union_atoms": union_atoms,
            "union_nitrogens": union_nitrogens,
            "ring_sizes": ring_sizes,
            "comp_indices": comp  # indices of rings in the candidate
        })
    
    # Now search among candidate fused systems for one that meets nucleobase-like criteria.
    candidate_found = None
    candidate_type = None
    # We require that the candidate system be the only aromatic heterocycle in the molecule.
    # In other words, if there is more than one disconnected fused system that contains N atoms, reject.
    if len(candidates) > 1:
        return False, "Molecule contains more than one independent aromatic system; expecting a single nucleobase-like core"
        
    # There should be exactly one connected component (or we choose the largest if many exist)
    candidate = candidates[0]
    
    # Check for pyrimidine-like (monocyclic, six-membered with at least 2 nitrogens)
    if candidate["ring_count"] == 1:
        # Retrieve the sole ring info.
        ring = rings[candidate["comp_indices"][0]]
        if ring["size"] == 6 and ring["n_count"] >= 2:
            candidate_type = "pyrimidine-like"
        else:
            return False, "The single aromatic ring does not meet criteria for a pyrimidine-like heterocycle (need six-membered with ≥2 nitrogens)"
    # Check for purine-like (fused system of exactly two rings, one six-membered and one five-membered,
    # sharing at least 2 atoms, and union has at least 3 nitrogens)
    elif candidate["ring_count"] == 2:
        # Identify the two rings.
        ring1 = rings[candidate["comp_indices"][0]]
        ring2 = rings[candidate["comp_indices"][1]]
        # Check sizes: one must be 6-membered and the other 5-membered.
        sizes = sorted([ring1["size"], ring2["size"]])
        if sizes != [5, 6]:
            return False, "The fused ring system is not composed of one 5-membered and one 6-membered ring (purine-like)"
        # They must share at least 2 atoms.
        shared = ring1["indices"].intersection(ring2["indices"])
        if len(shared) < 2:
            return False, "The two rings do not share at least 2 atoms; insufficient fusion for a purine-like system"
        # The union of the two rings must have at least 3 nitrogen atoms.
        if candidate["union_nitrogens"] < 3:
            return False, "The fused ring system does not have at least 3 nitrogen atoms required for a purine-like core"
        candidate_type = "purine-like"
    else:
        # More than two rings fused together (or less than 1) is unexpected for a canonical nucleobase analogue.
        return False, f"Candidate fused aromatic system contains {candidate['ring_count']} rings; expected 1 (pyrimidine-like) or 2 (purine-like)"
    
    # Finally, check that the candidate nucleobase core represents a significant fraction of the molecule’s heavy atoms.
    total_heavy = mol.GetNumHeavyAtoms()
    candidate_fraction = len(candidate["union_atoms"]) / total_heavy
    # We relax the fraction threshold to 0.35 rather than 0.5 to accommodate known nucleobase analogues with small substituents.
    if candidate_fraction < 0.35:
        return False, f"Candidate heterocycle accounts for only {candidate_fraction:.2f} fraction of heavy atoms (<0.35), suggesting too many extra substituents"
    
    return True, f"Molecule fits refined heuristic criteria as a nucleobase analogue ({candidate_type}); candidate core covers {candidate_fraction:.2f} of heavy atoms."

# Example usage:
if __name__ == "__main__":
    # You can test a few examples
    tests = {
        "N(4)-acetylcytosine": "O=C1NC=CC(=N1)NC(=O)C",
        "5-bromouracil": "Brc1c[nH]c(=O)[nH]c1=O",
        "8-azaadenine (false negative previously)": "Nc1ncnc2[nH]nnc12",
        "flucytosine (false positive previously)": "Nc1nc(=O)[nH]cc1F"
    }
    for name, smi in tests.items():
        res, reason = is_nucleobase_analogue(smi)
        print(f"{name}: {res}, {reason}")
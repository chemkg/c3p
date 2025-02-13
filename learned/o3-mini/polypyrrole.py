"""
Classifies: CHEBI:38077 polypyrrole
"""
#!/usr/bin/env python
"""
Classifies: A compound composed of two or more pyrrole units.
A pyrrole unit is defined here as a five‐membered aromatic ring 
that contains exactly one nitrogen atom and four carbon atoms.
To avoid false positives from isolated pyrrole fragments appended to a larger scaffold,
this version requires that at least two pyrrole matches are connected in the molecule.
Furthermore, we compute the fraction of heavy atoms (non‐hydrogens) that belong to the 
union of those connected pyrrole units. If exactly two connected rings are found we require
a lower fraction (>=20%), but if three or more connected rings are found we demand that they
contain at least 40% of the molecule’s heavy atoms.
Note that this rule‐based method is imperfect – some known polypyrroles may be missed and
some false positives may still pass.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole – i.e. it contains two or more pyrrole units in a connected network.
    A pyrrole unit is a five-membered aromatic ring containing one nitrogen atom and four carbon atoms.
    
    Our algorithm uses the following steps:
      1. Parse and sanitize the molecule.
      2. Count heavy atoms (non‐hydrogen atoms).
      3. Find all substructures that match a pyrrole SMARTS.
      4. If fewer than 2 pyrrole matches are found, return False.
      5. For each match, we "expand" its atom indices by including all direct neighbors.
      6. We then cluster pyrrole matches into connected (or near‐connected) groups.
         Two matches are considered connected if the expanded neighborhood of one overlaps 
         with the actual atom set of the other.
      7. For each connected cluster we compute (a) the number ("n_ring") of pyrrole units 
         and (b) the fraction (f) = (# unique heavy atoms in all pyrrole rings in the cluster)/(# heavy atoms in the molecule)
      8. Finally, we require that at least one cluster has:
           - n_ring >= 2 and,
           - if n_ring == 2: f >= 0.20, else (n_ring >= 3): f >= 0.40.
           
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule qualifies as a (connected) polypyrrole, False otherwise.
       str: Detailed reason for the classification decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule (to assign aromaticity etc.)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Count heavy (non-hydrogen) atoms in the molecule.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms found in molecule"
    total_heavy = len(heavy_atoms)
    
    # Define the SMARTS for a pyrrole unit. This pattern matches a 5-membered aromatic ring
    # with one nitrogen ([n;a,r]) and four aromatic carbons ([c;a,r]).
    pyrrole_smarts = "[n;a,r]1[c;a,r][c;a,r][c;a,r][c;a,r]1"
    pyrrole_query = Chem.MolFromSmarts(pyrrole_smarts)
    if pyrrole_query is None:
        return None, None  # in case SMARTS fails to compile
    
    # Find all matches of the pyrrole substructure.
    # Each match is a tuple of atom indices.
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_query, uniquify=True)
    num_matches = len(pyrrole_matches)
    if num_matches < 2:
        return False, f"Found {num_matches} pyrrole unit(s); at least 2 are required."
    
    # For each pyrrole match, construct an "expanded" set:
    # the atoms in the match plus any direct neighbor (bond-connected) atoms.
    expanded_matches = []
    for match in pyrrole_matches:
        atoms_set = set(match)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                atoms_set.add(nbr.GetIdx())
        expanded_matches.append(atoms_set)
    
    # Cluster the pyrrole matches into connected groups.
    # We use a simple union-find logic over the indices of pyrrole_matches.
    parent = list(range(num_matches))
    
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    
    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri
    
    # Two pyrrole matches are considered connected if the expanded set of one 
    # overlaps the original atom set of the other.
    for i in range(num_matches):
        for j in range(i+1, num_matches):
            set_i = set(pyrrole_matches[i])
            set_j = set(pyrrole_matches[j])
            # Check if the expanded neighbourhood of match i touches match j or vice versa.
            if (expanded_matches[i] & set_j) or (expanded_matches[j] & set_i):
                union(i, j)
    
    # Group the matches by cluster root.
    clusters = {}
    for i in range(num_matches):
        root = find(i)
        clusters.setdefault(root, []).append(i)
    
    # For each cluster, compute:
    #  a) the number of pyrrole units (n_ring),
    #  b) the union of atom indices (only those in the actual pyrrole rings, not the expansion).
    cluster_info = []
    for clust in clusters.values():
        n_ring = len(clust)
        union_atoms = set()
        for idx in clust:
            union_atoms.update(pyrrole_matches[idx])
        cluster_info.append((n_ring, union_atoms))
    
    # Choose the largest cluster (by number of rings) that is a candidate.
    best = None
    for n_ring, atoms_set in cluster_info:
        fraction = len(atoms_set) / total_heavy
        if best is None or n_ring > best[0]:
            best = (n_ring, atoms_set, fraction)
    
    # If no cluster found:
    if best is None:
        return False, "No connected pyrrole clusters found."
    
    n_ring, atoms_set, fraction = best
    # Set threshold based on the size of the connected cluster.
    if n_ring == 2:
        threshold = 0.20
    else:
        threshold = 0.40
    
    reason = (f"Connected cluster has {n_ring} pyrrole unit(s) covering {fraction:.2f} "
              f"of heavy atoms (threshold {threshold}); ")
    
    if n_ring < 2:
        return False, reason + "At least 2 connected pyrrole units are required."
    if fraction < threshold:
        return False, reason + "Molecule is not sufficiently composed of pyrrole units."
    
    return True, reason + "Qualifies as a polypyrrole."

# Example usage:
if __name__ == "__main__":
    # Here are a few examples.
    test_cases = [
        # A true positive example (porphyrin-like: several connected pyrrole units).
        ("Cc1c2Cc3[nH]c(Cc4[nH]c(Cc5[nH]c(Cc([nH]2)c1CCC(O)=O)c(C)c5CCC(O)=O)c(C)c4CCC(O)=O)c(C)c3CCC(O)=O",
         "coproporphyrinogen III"),
        # A dipyrromethane example.
        ("C(c1ccc[nH]1)c1ccc[nH]1", "dipyrromethane"),
        # Example that previously was a false negative: (1R)-primary fluorescent chlorophyll catabolite.
        ("C=1(C2=C(NC1CC=3NC(C(=O)[H])=C(C)C3CC)/C(=C/4\\N=C(C[C@@]5(C(C)=C(C=C)C(=O)N5)[H])[C@@H](C)[C@@H]4CCC(=O)O)/[C@@H](C(=O)OC)C2=O)",
         "(1R)-primary fluorescent chlorophyll catabolite"),
        # Example that was a false positive in the previous version.
        ("ClC1=C(Cl)C=C2C(C=3C(=C(C(=O)OC)NC3)C=4C5=C(C=CC(=C5)Cl)NC4)=CNC2=C1", "Lynamicin B")
    ]
    
    for smi, name in test_cases:
        result, reason = is_polypyrrole(smi)
        print(f"Name: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")
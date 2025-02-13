"""
Classifies: CHEBI:25676 oligopeptide
"""
#!/usr/bin/env python
"""
Classifies: A peptide containing a relatively small number of amino acids (oligopeptide)

This improved approach first scans for candidate peptide bonds by analyzing every single, non‐ring bond 
between a carbon and a nitrogen such that (i) the carbon (carbonyl candidate) has a double‐bonded O neighbor 
and (ii) the nitrogen (amide candidate) is attached to at least one carbon (an alpha carbon). 

All candidate bonds (stored as (carbon_index, nitrogen_index)) are then assembled into connected components 
(via an undirected graph) representing contiguous peptide chains. For a valid linear peptide, we require that 
(i) there is a free (terminal) carbonyl “end” (a carbon that appears only as the carbonyl in candidate bonds) 
and (ii) a free terminal amine (a nitrogen appearing only as such). 

We also check that the substructure in the candidate chain has a minimal molecular weight to somewhat reject 
small, chance amide bonds. Finally, the chain length is taken as the number of candidate bonds plus one, and we 
accept the chain only if its residue count is between 2 and 10.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    
    The algorithm:
      (1) Parse the molecule.
      (2) Iterate over every single bond (outside rings) to look for a bond between a carbon and a nitrogen
          that is consistent with a peptide bond, i.e.:
              - The carbon (carbonyl candidate) must have at least one neighbor oxygen connected via a double bond.
              - The nitrogen (amide candidate) must have at least one neighboring carbon (an alpha carbon) not being the carbonyl.
      (3) Each candidate bond is stored as (carbonyl_index, nitrogen_index).
      (4) An undirected connectivity graph is built from all atoms that occur in these candidate bonds.
      (5) The graph is partitioned into connected components. In a valid peptide chain, the candidate bonds occur
          in series so that (a) one carbon appears only as the carbonyl (N‐terminal “end” of a residue) and (b)
          one nitrogen appears only as the amide (C–terminal “end”). These terminal ideas are used to filter out
          spurious matches.
      (6) In any given connected component the chain length is defined as the number of candidate bonds plus one.
      (7) The substructure is extracted and its ExactMolWt is computed – a minimum weight threshold is applied 
          to further reduce false positives.
          
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is an oligopeptide (2–10 amino acids in a contiguous linear chain), False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    candidate_bonds = []  # list of tuples: (carbonyl_atom_idx, amide_nitrogen_idx)

    # Iterate over every bond in the molecule.
    for bond in mol.GetBonds():
        # Only consider single bonds that are not in a ring.
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        # Check for a bond between carbon and nitrogen.
        # We want the carbon to be a carbonyl carbon (i.e. has at least one double‐bonded oxygen neighbor)
        # and the other atom to be nitrogen that is connected to an alpha-carbon.
        carbon = None
        nitrogen = None
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            carbon = a1
            nitrogen = a2
        elif a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            carbon = a2
            nitrogen = a1
        else:
            continue

        # For the carbon candidate, check if it has an oxygen neighbor by a double bond.
        has_carbonyl = False
        for nb in carbon.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                bond_CO = mol.GetBondBetweenAtoms(carbon.GetIdx(), nb.GetIdx())
                if bond_CO and bond_CO.GetBondType() == Chem.BondType.DOUBLE:
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            continue

        # For the nitrogen candidate, check if it has at least one neighboring carbon
        # (this is our rough check that it is attached to an α-carbon).
        has_alpha = False
        for nb in nitrogen.GetNeighbors():
            if nb.GetIdx() == carbon.GetIdx():
                continue
            if nb.GetAtomicNum() == 6:
                has_alpha = True
                break
        if not has_alpha:
            continue

        # If both criteria are met, record the candidate peptide bond.
        candidate_bonds.append((carbon.GetIdx(), nitrogen.GetIdx()))

    if not candidate_bonds:
        return False, "No candidate peptide bonds detected; not a peptide."

    # Build an undirected connectivity graph for all atoms in candidate bonds.
    graph = {}
    for bond_tuple in candidate_bonds:
        for atom_idx in bond_tuple:
            if atom_idx not in graph:
                graph[atom_idx] = set()
        # Link the two atoms.
        graph[bond_tuple[0]].add(bond_tuple[1])
        graph[bond_tuple[1]].add(bond_tuple[0])

    # Partition the graph into connected components.
    visited = set()
    components = []
    for node in graph:
        if node in visited:
            continue
        comp = set()
        stack = [node]
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp.add(cur)
            for neighbor in graph[cur]:
                if neighbor not in comp:
                    stack.append(neighbor)
        components.append(comp)
        visited |= comp

    # For each component, count how many candidate bonds (within that component) appear and also verify some terminal conditions.
    max_chain_length = 0
    for comp in components:
        # Count how many candidate bonds are fully contained in this component.
        chain_bonds = [bond for bond in candidate_bonds if bond[0] in comp and bond[1] in comp]
        chain_length = len(chain_bonds) + 1  # residue count = peptide bonds + 1

        # Terminal conditions:
        # * Terminal carbon atoms: those that act as the carbonyl in candidate bonds and never as an amide nitrogen.
        # * Terminal nitrogen atoms: those that act as the amide and never as a carbonyl.
        carbons = [bond[0] for bond in candidate_bonds if bond[0] in comp]
        nitrogens = [bond[1] for bond in candidate_bonds if bond[1] in comp]
        terminal_carbons = set(carbons) - set(nitrogens)
        terminal_nitrogens = set(nitrogens) - set(carbons)
        if not terminal_carbons or not terminal_nitrogens:
            # The component does not display the expected free termini.
            continue

        # Extract the substructure corresponding to the component.
        try:
            submol = Chem.PathToSubmol(mol, list(comp))
        except Exception:
            # Could not extract the substructure; skip this component.
            continue
        mw = Descriptors.ExactMolWt(submol)
        # Assume that even a dipeptide should have a modest weight (e.g., more than 150 Da).
        if mw < 150:
            continue

        if chain_length > max_chain_length:
            max_chain_length = chain_length

    if max_chain_length == 0:
        return False, "No valid contiguous peptide chain detected."
    if max_chain_length < 2:
        return False, "Only a single residue detected; not a peptide."
    if max_chain_length > 10:
        return False, f"Peptide has {max_chain_length} amino acids, which is too many to be considered an oligopeptide."

    return True, f"Peptide detected with {max_chain_length} amino acids; classified as an oligopeptide."

# Uncomment below to test an example.
# test_smiles = "O=C(N[C@@H](CC1=CC=CC=C1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CCCN=C(N)N"  # Arg-Arg-Phe
# print(is_oligopeptide(test_smiles))
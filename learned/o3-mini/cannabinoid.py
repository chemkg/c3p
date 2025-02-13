"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid
Cannabinoids are a diverse group of pharmacologically active secondary metabolites 
that occur in the Cannabis plant and are also produced endogenously in humans and animals.
They are characterized by the presence of oxygen – either as part of a heterocyclic ring,
or in the form of various functional groups such as aromatic hydroxyls or ethanolamine/glycerol substructures.
This heuristic approach checks for one or more of these features.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    
    This heuristic algorithm requires that:
      1. The SMILES can be parsed.
      2. The molecule has at least one oxygen atom.
      3. The molecule displays at least one cannabinoid-like motif:
           - A ring containing an oxygen (as in benzopyran or related heterocycles),
           - OR an aromatic hydroxyl group (as in resorcinol moieties),
           - OR an ethanolamine/glycerol fragment (observed in endocannabinoids).
      4. The molecular weight is not unrealistically low (here we expect >200 Da).

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is likely a cannabinoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of oxygen atoms
    oxy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxy_atoms:
        return False, "Molecule has no oxygen atoms and thus cannot be a cannabinoid"

    # Check molecular weight (a crude filter; many cannabinoids are >200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a typical cannabinoid"

    found_cannabinoid_motif = False
    reason_parts = []

    # Check 1: Look for a heterocyclic ring containing oxygen.
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        # If the ring has at least one oxygen and at least one carbon, consider it.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetAtomicNum() == 8 for atom in ring_atoms) and any(atom.GetAtomicNum() == 6 for atom in ring_atoms):
            found_cannabinoid_motif = True
            reason_parts.append("found oxygen-containing heterocyclic ring")
            break

    # Check 2: Look for an aromatic hydroxyl motif.
    if not found_cannabinoid_motif:
        # This SMARTS finds an aromatic carbon bound directly to a hydroxyl group.
        arom_hydroxyl_smarts = Chem.MolFromSmarts("c[OX2H]")
        if mol.HasSubstructMatch(arom_hydroxyl_smarts):
            found_cannabinoid_motif = True
            reason_parts.append("found aromatic hydroxyl group (resorcinol-like motif)")

    # Check 3: Look for an ethanolamine fragment often seen in endocannabinoids.
    if not found_cannabinoid_motif:
        ethanolamine_smarts = Chem.MolFromSmarts("NCCO")
        if mol.HasSubstructMatch(ethanolamine_smarts):
            found_cannabinoid_motif = True
            reason_parts.append("found ethanolamine/glycerol fragment")

    if not found_cannabinoid_motif:
        return False, "Molecule does not display any cannabinoid-like oxygen motifs (no heterocyclic oxygen, aromatic hydroxyl, or ethanolamine fragment found)"
    
    # If any motif was found, we consider it a cannabinoid.
    combined_reason = "; ".join(reason_parts) + f"; molecular weight = {mol_wt:.1f} Da"
    return True, combined_reason
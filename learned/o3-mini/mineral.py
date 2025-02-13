"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: minerals

A mineral is normally an inorganic (often ionic) substance. Because properties such as 
crystallinity or geological history are not encoded in a SMILES string, the classifier uses 
structural proxies.

Heuristics used include:
  - Reject any species that contains an “organic” fragment (here defined as having at least 3 carbons 
    and 3 hydrogens in one connected fragment).
  - If there are multiple disconnected fragments then at least one fragment must contain a metal.
       Furthermore, if all “non–metal” fragments consist of a single atom that is a halide (Cl, Br, F or I)
       and the metal fragments are comprised of transition metals, then we reject the species (e.g. ZnCl2).
  - For a single–fragment molecule, we require that it contain a metal. We also reject if the molecule 
    contains a peroxo bond ([O]-[O]), which often flags peroxo or other reactive motifs rather than a 
    stable salt.
 
Note: These rules are necessarily heuristic and will fail on borderline cases.
"""

from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string using
    several heuristic checks.
    
    The algorithm works as follows:
      1. Parse the SMILES string into a molecule and disconnect into fragments.
      2. Reject if any fragment appears organic (>=3 C and >=3 H atoms after adding H's).
      3. For molecules with multiple disconnected fragments:
           - At least one fragment must contain a metal.
           - If all non‐metal fragments are single atoms and that metal (or metals) belong to classical 
             transition metals, then (for example in simple metal–halide salts like ZnCl2) the species 
             is judged non–mineral.
      4. For a single connected molecule:
           - If a peroxo bond ([O]-[O]) is found, reject.
           - Otherwise, require the molecule to contain at least one metal.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if classified as a mineral, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get disconnected fragments as separate molecules.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    
    # A helper function: Determine if a fragment is clearly organic.
    def is_fragment_organic(frag):
        # Add hydrogens for reliable atom counting.
        frag_with_H = Chem.AddHs(frag)
        count_C = sum(1 for atom in frag_with_H.GetAtoms() if atom.GetSymbol() == "C")
        count_H = sum(1 for atom in frag_with_H.GetAtoms() if atom.GetSymbol() == "H")
        return (count_C >= 3 and count_H >= 3)
    
    for frag in fragments:
        if is_fragment_organic(frag):
            return False, ("Contains organic fragment(s) (at least 3 C and 3 H atoms in a fragment), "
                           "indicating an organic salt rather than a typical inorganic mineral.")
    
    # Define sets of metals and halide counter-ions.
    metal_set = {
        # Alkali and alkaline earth:
        "Li","Na","K","Rb","Cs","Fr","Be","Mg","Ca","Sr","Ba","Ra",
        # Transition metals and some post-transition:
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Pd","Ag","Cd","Pt","Au","Hg",
        # Other common metals in minerals:
        "Al","Sb","Pb","Bi","Sn","In","Tl",
        # Lanthanides & actinides:
        "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
        "Ac","Th","Pa","U","Np","Pu"
    }
    transition_metals = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Pd","Ag","Cd","Pt","Au","Hg"}
    halide_set = {"F", "Cl", "Br", "I"}
    
    # For multi-fragment molecules:
    if len(fragments) > 1:
        # Check if at least one fragment contains a metal.
        has_metal = any(any(atom.GetSymbol() in metal_set for atom in frag.GetAtoms()) for frag in fragments)
        if not has_metal:
            return False, "Multiple fragments but none contain metal, not typical of inorganic minerals."
        
        # Check if the non-metal fragments are very simple (i.e. single atoms that are halides)
        non_metal_frags = [frag for frag in fragments if not any(atom.GetSymbol() in metal_set for atom in frag.GetAtoms())]
        if non_metal_frags and all(frag.GetNumAtoms() == 1 and frag.GetAtomWithIdx(0).GetSymbol() in halide_set
                                     for frag in non_metal_frags):
            # In this case, if the metal fragments are composed solely of transition metals,
            # we judge the salt as too simple (e.g. ZnCl2) and not assign it as a mineral.
            metal_frags = [frag for frag in fragments if any(atom.GetSymbol() in metal_set for atom in frag.GetAtoms())]
            if metal_frags and all(
                any(atom.GetSymbol() in transition_metals for atom in frag.GetAtoms())
                for frag in metal_frags
            ):
                return False, "Simple metal-halide salt of a transition metal, not typically classified as a mineral."
        
        return True, "Contains multiple disconnected fragments with metal, consistent with ionic compounds found in minerals."
    
    # For single–fragment molecules:
    else:
        # If an O–O bond is found (common in peroxo species), then reject.
        peroxo_pat = Chem.MolFromSmarts("[O]-[O]")
        if mol.HasSubstructMatch(peroxo_pat):
            return False, "Contains peroxo (O-O) bonds not typical of stable mineral phases."
        
        # For a single fragment, require a metal to be present.
        if any(atom.GetSymbol() in metal_set for atom in mol.GetAtoms()):
            return True, "Single fragment containing metal, consistent with many inorganic minerals."
        else:
            return False, "Single fragment with no metal, not typical of inorganic minerals."

# Example usage (for testing):
if __name__ == '__main__':
    test_compounds = [
        # True positives:
        ("[Na+].[Na+].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2", "disodium tetraborate octahydrate"),
        ("[Fe++].[S-][S-]", "pyrite"),
        ("[OH-].[OH-].[Ca++]", "calcium dihydroxide"),
        ("[O-]S([O-])(=O)=O.[Ba+2]", "barium sulfate"),
        ("[Mg++].[O-]S([O-])(=O)=O", "magnesium sulfate"),
        ("[Ca++].[H]OP([O-])([O-])=O", "calcium hydrogenphosphate"),
        ("[Pd-2](Cl)(Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+]", "Potassium hexachloropalladate(IV)"),
        ("[Ni]=S=[Ni]=S=[Ni]", "heazlewoodite"),
        ("[Fe++].[O-]C([O-])=O", "ferrous carbonate"),
        ("[OH-].[OH-].[OH-].[OH-].O=[Si]([O-])O[Si](=O)[O-].[Al+3].[Al+3]", "kaolinite"),
        ("O1B(O[B-]2(OB(O[B-]1(O[H])O2)O[H])O[H])O[H].[Na+].[Na+].O.O.O.O.O.O.O.O.O.O", "disodium tetraborate decahydrate"),
        ("[Ca++].[O-]S([O-])(=O)=O", "calcium sulfate"),
        ("[Cs+].[H]C([O-])=O", "caesium formate"),
        ("[Cl-].[Cl-].[Ca++]", "calcium dichloride"),
        ("[Ba++].CC([O-])=O.CC([O-])=O", "barium acetate"),
        ("O.O.Cl[Cu]Cl", "copper(II) chloride dihydrate"),
        ("[Sb](F)(F)(F)(F)F", "antimony pentafluoride"),
        ("[Na+].[Na+].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2", "disodium tetraborate"),
        ("[S--].[S--].[S--].[S--].[Fe++].[Fe+3].[Fe+3]", "greigite"),
        ("[Ba++].[O-]C([O-])=O", "barium carbonate"),
        ("O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O", "aluminium sulfate octadecahydrate"),
        ("[Zn++].[O-][N+]([O-])=O.[O-][N+]([O-])=O", "zinc nitrate"),
        ("[Ca+2].C(=O)([O-])[O-]", "calcium carbonate"),
        ("O[Ca]", "calcium monohydroxide"),
        ("Cl[Cu]Cl.O.O.O.O.O", "copper(II) chloride pentahydrate"),
        ("[Fe+3].[O-]P([O-])(=O)[O-]", "iron(3+) phosphate"),
        # The following should be rejected (false negatives in previous version):
        ("[Cl-].[Cs+]", "caesium chloride"),
        ("[Cl-].[K+]", "potassium chloride"),
        # Organic salt examples (should be rejected):
        ("[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O", "magnesium distearate"),
        ("[Mg++].CCC([O-])=O.CCC([O-])=O", "magnesium dipropionate"),
        # Some cases that should be rejected because they show peroxo/other motifs:
        ("O=[I+]=O", "dioxidoiodine"),
        ("[O-]OS([O-])(=O)=O", "peroxysulfate"),
        ("F[N+](F)(F)F", "tetrafluoroammonium"),
        ("[O-]c1c(Cl)c(Cl)c([O-])c(Cl)c1Cl", "halogenated benzene-olate"),
        ("[Cl-].[Cl-].[Zn++]", "zinc dichloride")
    ]
    
    for smi, name in test_compounds:
        result, reason = is_mineral(smi)
        print("SMILES:", smi)
        print("NAME:", name)
        print("is_mineral:", result)
        print("Reason:", reason)
        print("------")
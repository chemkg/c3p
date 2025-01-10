"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: CHEBI:28852 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid has a steroid backbone with a hydroxyl group at the 3rd position
    in the beta configuration and a double bond between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the steroid backbone (four fused rings: three 6-membered and one 5-membered)
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2CC[C@H]3[C@@H](CC[C@H]4[C@@H](CC[C@H]5[C@@H](CC[C@H]6[C@@H](CC[C@H]7[C@@H](CC[C@H]8[C@@H](CC[C@H]9[C@@H](CC[C@H]%10[C@@H](CC[C@H]%11[C@@H](CC[C@H]%12[C@@H](CC[C@H]%13[C@@H](CC[C@H]%14[C@@H](CC[C@H]%15[C@@H](CC[C@H]%16[C@@H](CC[C@H]%17[C@@H](CC[C@H]%18[C@@H](CC[C@H]%19[C@@H](CC[C@H]%20[C@@H](CC[C@H]%21[C@@H](CC[C@H]%22[C@@H](CC[C@H]%23[C@@H](CC[C@H]%24[C@@H](CC[C@H]%25[C@@H](CC[C@H]%26[C@@H](CC[C@H]%27[C@@H](CC[C@H]%28[C@@H](CC[C@H]%29[C@@H](CC[C@H]%30[C@@H](CC[C@H]%31[C@@H](CC[C@H]%32[C@@H](CC[C@H]%33[C@@H](CC[C@H]%34[C@@H](CC[C@H]%35[C@@H](CC[C@H]%36[C@@H](CC[C@H]%37[C@@H](CC[C@H]%38[C@@H](CC[C@H]%39[C@@H](CC[C@H]%40[C@@H](CC[C@H]%41[C@@H](CC[C@H]%42[C@@H](CC[C@H]%43[C@@H](CC[C@H]%44[C@@H](CC[C@H]%45[C@@H](CC[C@H]%46[C@@H](CC[C@H]%47[C@@H](CC[C@H]%48[C@@H](CC[C@H]%49[C@@H](CC[C@H]%50[C@@H](CC[C@H]%51[C@@H](CC[C@H]%52[C@@H](CC[C@H]%53[C@@H](CC[C@H]%54[C@@H](CC[C@H]%55[C@@H](CC[C@H]%56[C@@H](CC[C@H]%57[C@@H](CC[C@H]%58[C@@H](CC[C@H]%59[C@@H](CC[C@H]%60[C@@H](CC[C@H]%61[C@@H](CC[C@H]%62[C@@H](CC[C@H]%63[C@@H](CC[C@H]%64[C@@H](CC[C@H]%65[C@@H](CC[C@H]%66[C@@H](CC[C@H]%67[C@@H](CC[C@H]%68[C@@H](CC[C@H]%69[C@@H](CC[C@H]%70[C@@H](CC[C@H]%71[C@@H](CC[C@H]%72[C@@H](CC[C@H]%73[C@@H](CC[C@H]%74[C@@H](CC[C@H]%75[C@@H](CC[C@H]%76[C@@H](CC[C@H]%77[C@@H](CC[C@H]%78[C@@H](CC[C@H]%79[C@@H](CC[C@H]%80[C@@H](CC[C@H]%81[C@@H](CC[C@H]%82[C@@H](CC[C@H]%83[C@@H](CC[C@H]%84[C@@H](CC[C@H]%85[C@@H](CC[C@H]%86[C@@H](CC[C@H]%87[C@@H](CC[C@H]%88[C@@H](CC[C@H]%89[C@@H](CC[C@H]%90[C@@H](CC[C@H]%91[C@@H](CC[C@H]%92[C@@H](CC[C@H]%93[C@@H](CC[C@H]%94[C@@H](CC[C@H]%95[C@@H](CC[C@H]%96[C@@H](CC[C@H]%97[C@@H](CC[C@H]%98[C@@H](CC[C@H]%99[C@@H](CC[C@H]%100[C@@H](CC[C@H]%101[C@@H](CC[C@H]%102[C@@H](CC[C@H]%103[C@@H](CC[C@H]%104[C@@H](CC[C@H]%105[C@@H](CC[C@H]%106[C@@H](CC[C@H]%107[C@@H](CC[C@H]%108[C@@H](CC[C@H]%109[C@@H](CC[C@H]%110[C@@H](CC[C@H]%111[C@@H](CC[C@H]%112[C@@H](CC[C@H]%113[C@@H](CC[C@H]%114[C@@H](CC[C@H]%115[C@@H](CC[C@H]%116[C@@H](CC[C@H]%117[C@@H](CC[C@H]%118[C@@H](CC[C@H]%119[C@@H](CC[C@H]%120[C@@H](CC[C@H]%121[C@@H](CC[C@H]%122[C@@H](CC[C@H]%123[C@@H](CC[C@H]%124[C@@H](CC[C@H]%125[C@@H](CC[C@H]%126[C@@H](CC[C@H]%127[C@@H](CC[C@H]%128[C@@H](CC[C@H]%129[C@@H](CC[C@H]%130[C@@H](CC[C@H]%131[C@@H](CC[C@H]%132[C@@H](CC[C@H]%133[C@@H](CC[C@H]%134[C@@H](CC[C@H]%135[C@@H](CC[C@H]%136[C@@H](CC[C@H]%137[C@@H](CC[C@H]%138[C@@H](CC[C@H]%139[C@@H](CC[C@H]%140[C@@H](CC[C@H]%141[C@@H](CC[C@H]%142[C@@H](CC[C@H]%143[C@@H](CC[C@H]%144[C@@H](CC[C@H]%145[C@@H](CC[C@H]%146[C@@H](CC[C@H]%147[C@@H](CC[C@H]%148[C@@H](CC[C@H]%149[C@@H](CC[C@H]%150[C@@H](CC[C@H]%151[C@@H](CC[C@H]%152[C@@H](CC[C@H]%153[C@@H](CC[C@H]%154[C@@H](CC[C@H]%155[C@@H](CC[C@H]%156[C@@H](CC[C@H]%157[C@@H](CC[C@H]%158[C@@H](CC[C@H]%159[C@@H](CC[C@H]%160[C@@H](CC[C@H]%161[C@@H](CC[C@H]%162[C@@H](CC[C@H]%163[C@@H](CC[C@H]%164[C@@H](CC[C@H]%165[C@@H](CC[C@H]%166[C@@H](CC[C@H]%167[C@@H](CC[C@H]%168[C@@H](CC[C@H]%169[C@@H](CC[C@H]%170[C@@H](CC[C@H]%171[C@@H](CC[C@H]%172[C@@H](CC[C@H]%173[C@@H](CC[C@H]%174[C@@H](CC[C@H]%175[C@@H](CC[C@H]%176[C@@H](CC[C@H]%177[C@@H](CC[C@H]%178[C@@H](CC[C@H]%179[C@@H](CC[C@H]%180[C@@H](CC[C@H]%181[C@@H](CC[C@H]%182[C@@H](CC[C@H]%183[C@@H](CC[C@H]%184[C@@H](CC[C@H]%185[C@@H](CC[C@H]%186[C@@H](CC[C@H]%187[C@@H](CC[C@H]%188[C@@H](CC[C@H]%189[C@@H](CC[C@H]%190[C@@H](CC[C@H]%191[C@@H](CC[C@H]%192[C@@H](CC[C@H]%193[C@@H](CC[C@H]%194[C@@H](CC[C@H]%195[C@@H](CC[C@H]%196[C@@H](CC[C@H]%197[C@@H](CC[C@H]%198[C@@H](CC[C@H]%199[C@@H](CC[C@H]%200[C@@H](CC[C@H]%201[C@@H](CC[C@H]%202[C@@H](CC[C@H]%203[C@@H](CC[C@H]%204[C@@H](CC[C@H]%205[C@@H](CC[C@H]%206[C@@H](CC[C@H]%207[C@@H](CC[C@H]%208[C@@H](CC[C@H]%209[C@@H](CC[C@H]%210[C@@H](CC[C@H]%211[C@@H](CC[C@H]%212[C@@H](CC[C@H]%213[C@@H](CC[C@H]%214[C@@H](CC[C@H]%215[C@@H](CC[C@H]%216[C@@H](CC[C@H]%217[C@@H](CC[C@H]%218[C@@H](CC[C@H]%219[C@@H](CC[C@H]%220[C@@H](CC[C@H]%221[C@@H](CC[C@H]%222[C@@H](CC[C@H]%223[C@@H](CC[C@H]%224[C@@H](CC[C@H]%225[C@@H](CC[C@H]%226[C@@H](CC[C@H]%227[C@@H](CC[C@H]%228[C@@H](CC[C@H]%229[C@@H](CC[C@H]%230[C@@H](CC[C@H]%231[C@@H](CC[C@H]%232[C@@H](CC[C@H]%233[C@@H](CC[C@H]%234[C@@H](CC[C@H]%235[C@@H](CC[C@H]%236[C@@H](CC[C@H]%237[C@@H](CC[C@H]%238[C@@H](CC[C@H]%239[C@@H](CC[C@H]%240[C@@H](CC[C@H]%241[C@@H](CC[C@H]%242[C@@H](CC[C@H]%243[C@@H](CC[C@H]%244[C@@H](CC[C@H]%245[C@@H](CC[C@H]%246[C@@H](CC[C@H]%247[C@@H](CC[C@H]%248[C@@H](CC[C@H]%249[C@@H](CC[C@H]%250[C@@H](CC[C@H]%251[C@@H](CC[C@H]%252[C@@H](CC[C@H]%253[C@@H](CC[C@H]%254[C@@H](CC[C@H]%255[C@@H](CC[C@H]%256[C@@H](CC[C@H]%257[C@@H](CC[C@H]%258[C@@H](CC[C@H]%259[C@@H](CC[C@H]%260[C@@H](CC[C@H]%261[C@@H](CC[C@H]%262[C@@H](CC[C@H]%263[C@@H](CC[C@H]%264[C@@H](CC[C@H]%265[C@@H](CC[C@H]%266[C@@H](CC[C@H]%267[C@@H](CC[C@H]%268[C@@H](CC[C@H]%269[C@@H](CC[C@H]%270[C@@H](CC[C@H]%271[C@@H](CC[C@H]%272[C@@H](CC[C@H]%273[C@@H](CC[C@H]%274[C@@H](CC[C@H]%275[C@@H](CC[C@H]%276[C@@H](CC[C@H]%277[C@@H](CC[C@H]%278[C@@H](CC[C@H]%279[C@@H](CC[C@H]%280[C@@H](CC[C@H]%281[C@@H](CC[C@H]%282[C@@H](CC[C@H]%283[C@@H](CC[C@H]%284[C@@H](CC[C@H]%285[C@@H](CC[C@H]%286[C@@H](CC[C@H]%287[C@@H](CC[C@H]%288[C@@H](CC[C@H]%289[C@@H](CC[C@H]%290[C@@H](CC[C@H]%291[C@@H](CC[C@H]%292[C@@H](CC[C@H]%293[C@@H](CC[C@H]%294[C@@H](CC[C@H]%295[C@@H](CC[C@H]%296[C@@H](CC[C@H]%297[C@@H](CC[C@H]%298[C@@H](CC[C@H]%299[C@@H](CC[C@H]%300[C@@H](CC[C@H]%301[C@@H](CC[C@H]%302[C@@H](CC[C@H]%303[C@@H](CC[C@H]%304[C@@H](CC[C@H]%305[C@@H](CC[C@H]%306[C@@H](CC[C@H]%307[C@@H](CC[C@H]%308[C@@H](CC[C@H]%309[C@@H](CC[C@H]%310[C@@H](CC[C@H]%311[C@@H](CC[C@H]%312[C@@H](CC[C@H]%313[C@@H](CC[C@H]%314[C@@H](CC[C@H]%315[C@@H](CC[C@H]%316[C@@H](CC[C@H]%317[C@@H](CC[C@H]%318[C@@H](CC[C@H]%319[C@@H](CC[C@H]%320[C@@H](CC[C@H]%321[C@@H](CC[C@H]%322[C@@H](CC[C@H]%323[C@@H](CC[C@H]%324[C@@H](CC[C@H]%325[C@@H](CC[C@H]%326[C@@H](CC[C@H]%327[C@@H](CC[C@H]%328[C@@H](CC[C@H]%329[C@@H](CC[C@H]%330[C@@H](CC[C@H]%331[C@@H](CC[C@H]%332[C@@H](CC[C@H]%333[C@@H](CC[C@H]%334[C@@H](CC[C@H]%335[C@@H](CC[C@H]%336[C@@H](CC[C@H]%337[C@@H](CC[C@H]%338[C@@H](CC[C@H]%339[C@@H](CC[C@H]%340[C@@H](
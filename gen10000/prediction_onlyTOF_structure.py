
#相関係数をプロット  https://smart-hint.com/python/corr/?msclkid=57501408cd2611eca9940947c336ccb5


from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem, Descriptors
from mordred import descriptors, Calculator
import numpy as np
import pandas as pd 
from pandas import DataFrame 
from sklearn.preprocessing import StandardScaler
from sklearn import model_selection
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.plotting import parallel_coordinates
import math
import  _pickle as cPickle
import gzip, numpy, copy, math
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import rdMolDescriptors
from matplotlib import cm
from sklearn.ensemble import RandomForestClassifier
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import rdMolDescriptors
from rdkit.Avalon import pyAvalonTools
from pandas import DataFrame as df
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import KFold

##########################
def make_df(smiles, mol):
    
    APR = []     # Aromatic Ratio
    AROM = []    # Aromatic Ring Count   芳香環の数（Cもヘテロも含む）
    CROM = []    # C芳香環の数
    HAROM = []   # ヘテロ環の数
    Ar_sp3 = []  # 芳香族原子数　–　sp3炭素数
    # BI = []      # C芳香環の数を全ての芳香環の数で除したもの
    oriBI=[]

    N=[]
    n=[]
    O=[]
    o=[]
    S=[]
    s=[]

    def Calc_ARR(mh):
        m = Chem.RemoveHs(mh)
        num_bonds = m.GetNumBonds()
        num_aromatic_bonds = 0
        for bond in m.GetBonds():
            if bond.GetIsAromatic():
                num_aromatic_bonds += 1
        ARR = num_aromatic_bonds/num_bonds
        return ARR


    A=Calc_ARR(mol)
    APR.append(A)
    apr=DataFrame({'APR':APR})

    
    #NAR AROM
    def Calc_Carbo_Hetero_Aromatic(mh):
        m = Chem.RemoveHs(mh)
        ring_info = m.GetRingInfo()
        atoms_in_rings = ring_info.AtomRings()
        num_Caromatic_ring = 0
        num_Hetaromatic_ring = 0
        num_Aromatic_ring = 0
        for ring in atoms_in_rings:
            aromatic_atom_in_ring = 0
            heteroatom_in_ring = 0
            for atom_id in ring:
                atom = m.GetAtomWithIdx(atom_id)
                if atom.GetIsAromatic():
                    aromatic_atom_in_ring += 1
                if atom.GetSymbol() != 'C': ### 環内の原子が炭素かどうかをチェック
                    heteroatom_in_ring += 1
            if aromatic_atom_in_ring == len(ring):
                if heteroatom_in_ring == 0:
                    num_Caromatic_ring += 1
                else:
                    num_Hetaromatic_ring += 1
            num_Aromatic_ring = num_Caromatic_ring + num_Hetaromatic_ring

        oBI = num_Caromatic_ring/num_Aromatic_ring

        return num_Aromatic_ring, num_Caromatic_ring, num_Hetaromatic_ring, oBI

   
    A, C, H, oBI= Calc_Carbo_Hetero_Aromatic(mol)
    AROM.append(A)
    CROM.append(C)
    HAROM.append(H)
    oriBI.append(oBI)

    arom=DataFrame({'AROM':AROM})
    crom=DataFrame({'CROM':CROM})
    harom=DataFrame({'HAROM':HAROM})
    oribi=DataFrame({'oriBI':oriBI})


    def Calc_Ar_Alk_balance(mh):
        m = Chem.RemoveHs(mh)
        num_aromatic_carbon = len(m.GetAromaticAtoms())
        num_sp3_carbon = 0
        for atom in m.GetAtoms():
            if str(atom.GetHybridization()) == 'SP3' and atom.GetSymbol() == 'C':
                num_sp3_carbon += 1
        ar_alk_balance = num_aromatic_carbon - num_sp3_carbon
        return ar_alk_balance



    A=Calc_Ar_Alk_balance(mol)
    Ar_sp3.append(A)

    ar_sp3=DataFrame({'Ar_sp3':Ar_sp3})





    #smi_list=list(df2['smiles'])

    #print(m)

    N_count=n_count=O_count=o_count=s_count=S_count=0    
    A=smiles
    #print(A)
    for x in range(len(A)):
        if A[x]=='N':
            N_count+=1
        if A[x]=='n':
            n_count+=1
        if A[x]=='O':
            O_count+=1
        if A[x]=='o':
            o_count+=1
        if A[x]=='S':
            S_count+=1
        if A[x]=='i':
            S_count+=-1   # Siを除く
        if A[x]=='s':
            s_count+=1
        if A[x]=='e':
            s_count+=-1   # snを除く

    N.append(N_count)
    n.append(n_count)
    O.append(O_count)
    o.append(o_count)
    S.append(S_count)
    s.append(s_count)

    num_N=DataFrame({'num_N':N})
    num_n=DataFrame({'num_n':n})
    num_O=DataFrame({'num_O':O})
    num_o=DataFrame({'num_o':o})
    num_S=DataFrame({'num_S':S})
    num_s=DataFrame({'num_s':s})


    #for m in ['CCCCN', 'c1ccccc1', 'c1ccncc1','Nc1ccncc1']:



    A=[]
    B=[]
    C=[]
    D=[]
    E=[]
    F=[]
    G=[]
    H=[]
    I=[]





    A.append(Chem.Descriptors.ExactMolWt(mol))  # print(type(A)) = list
    a=DataFrame({'MolWt':A})
    

    B.append(Chem.Descriptors.NumHeteroatoms(mol))
    b=DataFrame({'NumHeteroatoms':B})


    C.append(Chem.Descriptors.NumHAcceptors(mol))
    c=DataFrame({'NumHAcceptors':C})

    D.append(Chem.Descriptors.NumHDonors(mol))
    d=DataFrame({'NumHDonors':D})

    E.append(Chem.Descriptors.MolLogP(mol))
    e=DataFrame({'MolLogP':E})

    F.append(Chem.Descriptors.NumRotatableBonds(mol))
    f=DataFrame({'NumRotatableBonds':F})

    G.append(Chem.Descriptors.FractionCSP3(mol))
    g=DataFrame({'FractionCSP3':G})

    H.append(Chem.Descriptors.TPSA(mol))
    h=DataFrame({'TPSA':H})

    I.append(Chem.Descriptors.RingCount(mol))
    i=DataFrame({'RingCount':I})




    #説明変数X作成
    df4= pd.concat([a,b,c,d,e,f,g,h,i,apr,arom,crom,harom,ar_sp3,oribi,num_N,num_n,num_O,num_o,num_S,num_s], axis=1)
    

    return(df4)
###########################



def make_rf(smiles):
    dfz = PandasTools.LoadSDF('../datasets/Hole_mobility_onlyTOF300.sdf', molColName='ROMol',embedProps=True) 

    APR = []     # Aromatic Ratio
    AROM = []    # Aromatic Ring Count   芳香環の数（Cもヘテロも含む）
    CROM = []    # C芳香環の数
    HAROM = []   # ヘテロ環の数
    Ar_sp3 = []  # 芳香族原子数　–　sp3炭素数
    # BI = []      # C芳香環の数を全ての芳香環の数で除したもの
    oriBI=[]

    N=[]
    n=[]
    O=[]
    o=[]
    S=[]
    s=[]

    ''' apr,arom,crom,harom,ar_sp3,oribi,num_N,num_n,num_O,num_o,num_S,num_s '''
    #APR
    def Calc_ARR(mh):
        m = Chem.RemoveHs(mh)
        num_bonds = m.GetNumBonds()
        num_aromatic_bonds = 0
        for bond in m.GetBonds():
            if bond.GetIsAromatic():
                num_aromatic_bonds += 1
        ARR = num_aromatic_bonds/num_bonds
        return ARR

    for m in dfz.ROMol:
        A=Calc_ARR(m)
        APR.append(A)

    apr=DataFrame({'APR':APR})


    #NAR AROM
    def Calc_Carbo_Hetero_Aromatic(mh):
        m = Chem.RemoveHs(mh)
        ring_info = m.GetRingInfo()
        atoms_in_rings = ring_info.AtomRings()
        num_Caromatic_ring = 0
        num_Hetaromatic_ring = 0
        num_Aromatic_ring = 0
        for ring in atoms_in_rings:
            aromatic_atom_in_ring = 0
            heteroatom_in_ring = 0
            for atom_id in ring:
                atom = m.GetAtomWithIdx(atom_id)
                if atom.GetIsAromatic():
                    aromatic_atom_in_ring += 1
                if atom.GetSymbol() != 'C': ### 環内の原子が炭素かどうかをチェック
                    heteroatom_in_ring += 1
            if aromatic_atom_in_ring == len(ring):
                if heteroatom_in_ring == 0:
                    num_Caromatic_ring += 1
                else:
                    num_Hetaromatic_ring += 1
            num_Aromatic_ring = num_Caromatic_ring + num_Hetaromatic_ring

        oBI = num_Caromatic_ring/num_Aromatic_ring

        return num_Aromatic_ring, num_Caromatic_ring, num_Hetaromatic_ring, oBI

    for m in dfz.ROMol:
        A, C, H, oBI= Calc_Carbo_Hetero_Aromatic(m)
        AROM.append(A)
        CROM.append(C)
        HAROM.append(H)
        oriBI.append(oBI)

    arom=DataFrame({'AROM':AROM})
    crom=DataFrame({'CROM':CROM})
    harom=DataFrame({'HAROM':HAROM})
    oribi=DataFrame({'oriBI':oriBI})


    def Calc_Ar_Alk_balance(mh):
        m = Chem.RemoveHs(mh)
        num_aromatic_carbon = len(m.GetAromaticAtoms())
        num_sp3_carbon = 0
        for atom in m.GetAtoms():
            if str(atom.GetHybridization()) == 'SP3' and atom.GetSymbol() == 'C':
                num_sp3_carbon += 1
        ar_alk_balance = num_aromatic_carbon - num_sp3_carbon
        return ar_alk_balance


    for m in dfz.ROMol:
        A=Calc_Ar_Alk_balance(m)
        Ar_sp3.append(A)

    ar_sp3=DataFrame({'Ar_sp3':Ar_sp3})




    # mobility data
    df2=dfz.dropna(axis=1)
    smi_list=list(df2['smiles'])


    for m in smi_list:
        #print(m)
        mol=Chem.MolFromSmiles(m)
        N_count=n_count=O_count=o_count=s_count=S_count=0    
        A=Chem.MolToSmiles(mol)
        #print(A)
        for x in range(len(A)):
            if A[x]=='N':
                N_count+=1
            if A[x]=='n':
                n_count+=1
            if A[x]=='O':
                O_count+=1
            if A[x]=='o':
                o_count+=1
            if A[x]=='S':
                S_count+=1
            if A[x]=='i':
                S_count+=-1   # Siを除く
            if A[x]=='s':
                s_count+=1
            if A[x]=='e':
                s_count+=-1   # snを除く

        N.append(N_count)
        n.append(n_count)
        O.append(O_count)
        o.append(o_count)
        S.append(S_count)
        s.append(s_count)

        num_N=DataFrame({'num_N':N})
        num_n=DataFrame({'num_n':n})
        num_O=DataFrame({'num_O':O})
        num_o=DataFrame({'num_o':o})
        num_S=DataFrame({'num_S':S})
        num_s=DataFrame({'num_s':s})



    A=[]
    B=[]
    C=[]
    D=[]
    E=[]
    F=[]
    G=[]
    H=[]
    I=[]




    for m in dfz.ROMol:
        A.append(Chem.Descriptors.ExactMolWt(m))
        a=DataFrame({'MolWt':A})
    for m in dfz.ROMol:
        B.append(Chem.Descriptors.NumHeteroatoms(m))
        b=DataFrame({'NumHeteroatoms':B})

    for m in dfz.ROMol:
        C.append(Chem.Descriptors.NumHAcceptors(m))
        c=DataFrame({'NumHAcceptors':C})
    for m in dfz.ROMol:
        D.append(Chem.Descriptors.NumHDonors(m))
        d=DataFrame({'NumHDonors':D})
    for m in dfz.ROMol:
        E.append(Chem.Descriptors.MolLogP(m))
        e=DataFrame({'MolLogP':E})
    for m in dfz.ROMol:
        F.append(Chem.Descriptors.NumRotatableBonds(m))
        f=DataFrame({'NumRotatableBonds':F})
    for m in dfz.ROMol:
        G.append(Chem.Descriptors.FractionCSP3(m))
        g=DataFrame({'FractionCSP3':G})
    for m in dfz.ROMol:
        H.append(Chem.Descriptors.TPSA(m))
        h=DataFrame({'TPSA':H})
    for m in dfz.ROMol:
        I.append(Chem.Descriptors.RingCount(m))
        i=DataFrame({'RingCount':I})



    #Smiles読み込み
    mols = Chem.SDMolSupplier('../datasets/Hole_mobility_onlyTOF300.sdf')
    #フィンガープリントFCFP
    fps1 = [ AllChem.GetMorganFingerprintAsBitVect(mol,radius=2,nBits=1024) for mol in mols]   # ,useFeatures=True
    fps2 = [ list(map(int,list(fps))) for fps in fps1]
    fps=pd.DataFrame(fps2)



    #説明変数X作成
    df4= pd.concat([a,b,c,d,e,f,g,h,i,apr,arom,crom,harom,ar_sp3,oribi,num_N,num_n,num_O,num_o,num_S,num_s], axis=1)

    
    
    #sdfから必要な情報の抽出
    # 新規構造情報
    z3=pd.DataFrame(df4,columns=['MolWt','NumHeteroatoms','NumHAcceptors','NumHDonors','MolLogP','NumRotatableBonds','FractionCSP3','TPSA',
                                 'RingCount', 'APR', 'AROM', 'CROM', 'HAROM', 'oriBI', 'Ar_sp3', 'num_N', 'num_n', 'num_O', 'num_o', 'num_S', 'num_s'])
    #説明変数X作成
    z3= pd.concat([z3,fps], axis=1)         #　!!!!!!! ここのコメントを解除するとECFP考慮する。 !!!!!!!

    z4=z3[[898, 'num_N', 875, 139, 896, 'num_o', 215, 698, 212, 561, 956, 829, 891, 63, 736, 675, 950, 'NumHDonors', 19, 64, 785, 799, 33, 1013, 136, 338, 929, 'HAROM', 257, 549, 566, 261, 364, 73, 940, 114, 119, 15, 45, 271, 'num_n', 233, 801, 1009, 393, 'CROM', 723, 315, 'FractionCSP3', 452, 305, 655, 'NumHeteroatoms', 746, 'NumHAcceptors', 'NumRotatableBonds', 'oriBI', 378, 1017, 204, 'RingCount', 'Ar_sp3', 'MolWt', 'APR', 'MolLogP', 790, 'TPSA', 'AROM', 'num_S']]

    X = np.array(z4, dtype = np.float32)
    
    
    
    #目的変数ｙ作成
    sdf = [ mol for mol in Chem.SDMolSupplier('../datasets/Hole_mobility_onlyTOF300.sdf')]
    def getResponse( mols, prop= "mobility" ):
        y = []
        for mol in mols:
            act = mol.GetProp( prop )
            y.append( act )
        return y
    y = getResponse(sdf)
    y = np.array(y, dtype = np.float32)
    y = np.log10(y)
    y=np.reshape(y,(-1))  # yをi次元に。
    y = pd.DataFrame(y)


    #ランダムフォレスト
    from sklearn.model_selection import train_test_split, cross_val_score
    from sklearn.metrics import mean_squared_error
    from sklearn.ensemble import RandomForestRegressor


    #ランダム定数=rsでデータを0.75:0.25に分割    
    rsnum=591598
    X_train, X_test, y_train, y_test = train_test_split(X, y,test_size=0.20, random_state=rsnum)


    # regressor
    regressor= RandomForestRegressor(n_estimators=18, max_depth=14, max_features=38, random_state=2201)
    regressor.fit(X_train, y_train.values.ravel())
    
    #######################################################################################################
    # smiles 読み込み
    mol_smi =  Chem.MolFromSmiles(smiles)   # smilesをmolに変換
    df_smi=make_df(smiles, mol_smi)
    



    #フィンガープリントFCFP
    bit = {}        # フィンガープリントの情報を入れる空の辞書
    mol=mol_smi
    fps1 = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024, bitInfo=bit)  # , useFeatures=True
    fps2 = list(map(int,list(fps1)))
    fps=pd.DataFrame(fps2).T


    # smiles の説明変数を取得
    z = pd.concat([df_smi,fps], axis=1)
    z_smi = z[[898, 'num_N', 875, 139, 896, 'num_o', 215, 698, 212, 561, 956, 829, 891, 63, 736, 675, 950, 'NumHDonors', 19, 64, 785, 799, 33, 1013, 136, 338, 929, 'HAROM', 257, 549, 566, 261, 364, 73, 940, 114, 119, 15, 45, 271, 'num_n', 233, 801, 1009, 393, 'CROM', 723, 315, 'FractionCSP3', 452, 305, 655, 'NumHeteroatoms', 746, 'NumHAcceptors', 'NumRotatableBonds', 'oriBI', 378, 1017, 204, 'RingCount', 'Ar_sp3', 'MolWt', 'APR', 'MolLogP', 790, 'TPSA', 'AROM', 'num_S']]

    # smiles の説明変数Xを取得
    X_smi = np.array(z_smi, dtype = np.float32)


    # smiles のHMを予測
    ythree_HM = regressor.predict(X_smi)
    predict_31SR = 10**ythree_HM
    # np.set_printoptions(formatter={'float': '{:.2e}'.format})
    predict=float(predict_31SR)    # predict_31SRはnumpy.ndarray
    
    return(predict)

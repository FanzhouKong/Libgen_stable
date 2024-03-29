
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from rdkit import Chem
from molmass import Formula
from tqdm import tqdm
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from toolsets.search import string_search
# from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from toolsets.API_gets import inchi_to_smiles
def get_clean_std_list(std_list_raw, smile_col = 'smiles_fetched', info_col = None):
    if smile_col not in std_list_raw.columns and info_col is not None:
        std_list = complete_smiles(std_list_raw, reference_col=info_col)
    else:
        std_list = std_list_raw.copy()
    std_list = complete_formula(std_list, 'smiles_fetched')
    std_list_unsalted =neutrilize_salt_df(salt_df_input=std_list)
    std_list_cleaned = std_list_cleanup(std_list_unsalted, 'name')
    return(std_list_cleaned)

def std_list_cleanup(std_list_nt, name_column = 'name'):
    std_list = { 'name': std_list_nt[name_column],
             # 'abb':std_list_nt['abbreviated'],
            'inchikey': std_list_nt['uncharged_inchikey'],
            'mix':std_list_nt['mix'],
                 'smiles':std_list_nt['uncharged_smiles'],
                 'formula':std_list_nt['uncharged_formula'],
                 'mono_mass':std_list_nt['monoisotopic_mass'],
                 'formal_charges':std_list_nt['uncharged_formal_charges'],
             # 'rt_plus':std_list_nt['rt_plus'],
             # 'rt_minus':std_list_nt['rt_minus']
                 }
    std_list = pd.DataFrame(std_list)
    return(std_list)
# def get_info_from_db(std_list, inchikey_column = 'inchikey'):
#     smiles = []
#     std_list.dropna(axis=0, subset='inchikey', inplace=True)
#     for index, row in tqdm(std_list.iterrows(), total=len(std_list)):
#         try:
#             smiles.append(pubchem_get(content= row[inchikey_column]))
#         except:
#             smiles.append(np.NAN)
#     std_list['smiles_fetched']=smiles
#     # return(std_list)
#     # std_list.dropna(axis = 1, subset=['smiles_fetched'], inplace = True)
#     formal_charges = []
#     formula_fetched=[]
#     for index, row in std_list.iterrows():
#         try:
#             mol = Chem.MolFromSmiles(row['smiles_fetched'])
#             formula_fetched.append(CalcMolFormula(mol))
#             formal_charges.append(Chem.GetFormalCharge(mol))
#         except:
#             formula_fetched.append(np.NAN)
#             formal_charges.append(np.NAN)
#
#     std_list['formal_charges']=formal_charges
#     std_list['formula_fetched']=formula_fetched
#     std_list.dropna(axis=0, subset='smiles_fetched', inplace=True)
#     std_list.reset_index(inplace=True, drop=True)
#     return (std_list)



def complete_adducts(data, smile_column = 'smiles', adducts = ['[M+H]+']):
    for adduct in adducts:
        precursors = []
        for index, row in data.iterrows():
            precursors.append(calculate_precursormz(row[smile_column], adduct))
        data[adduct]= precursors
    return(data)


def check_mol(smile):
    if type(smile) is Chem.rdchem.Mol:
        mol_ = smile
        return(mol_)
    else:
        try:
            mol_ = Chem.MolFromSmiles(smile)
            return(mol_)
        except:
            return(np.NAN)
def cal_formal_charge(smile):
    mol_ = check_mol(smile)
    return (Chem.GetFormalCharge(mol_))

def check_salt(smile, salts_smart):
    mol_ = check_mol(smile)
    salt_pattern = Chem.MolFromSmarts(salts_smart)
    matches = mol_.GetSubstructMatches(salt_pattern)
    return(len(matches))
def neutrilize_salt(smile, salts_smart, return_stripped = False):
    mol_ = check_mol(smile)
    from rdkit.Chem import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize
    un = rdMolStandardize.Uncharger()
    remover = SaltRemover.SaltRemover(defnData=salts_smart)
    res, deleted = remover.StripMolWithDeleted(mol_)
    mol_uncharged = un.uncharge(res)
    if return_stripped == True:
        return(Chem.MolToSmiles(mol_uncharged),Chem.MolToSmiles(res))
    else:
        return(Chem.MolToSmiles(mol_uncharged))


def neutrilize_salt_df(salt_df_input, salts_smart='[Na+,K+,Cl-,HCl,HBr]', smile_column = 'smiles_fetched'):
    salt_df = salt_df_input.copy()
    salt_df.columns= salt_df.columns.str.lower()
    uncharged_result = {}
    for head in ['stripped','uncharged']:
        for tail in ['smiles', 'formula', 'formal_charges']:
            uncharged_result[head+'_'+tail]=[]
    uncharged_result['monoisotopic_mass']=[]
    for index, row in salt_df.iterrows():
        uncharged_smile,res_smile  = neutrilize_salt(row[smile_column], salts_smart, return_stripped=True)
        res = Chem.MolFromSmiles(res_smile)
        mol_uncharged = Chem.MolFromSmiles(uncharged_smile)
        for tail in ['smiles', 'formula', 'formal_charges','monoisotopic_mass']:
            if tail =='smiles':
                uncharged_result['stripped'+'_'+tail].append(Chem.MolToSmiles(res))
                uncharged_result['uncharged'+'_'+tail].append(Chem.MolToSmiles(mol_uncharged))
            if tail =='formula':
                uncharged_result['stripped'+'_'+tail].append(CalcMolFormula(res))
                uncharged_result['uncharged'+'_'+tail].append(CalcMolFormula(mol_uncharged))
            if tail =='formal_charges':
                uncharged_result['stripped'+'_'+tail].append(cal_formal_charge(res))
                uncharged_result['uncharged'+'_'+tail].append(cal_formal_charge(mol_uncharged))
            if tail == 'monoisotopic_mass':
                uncharged_result[tail].append(ExactMolWt(mol_uncharged))
    for i in uncharged_result.keys():
        salt_df[i]=uncharged_result[i]
    salt_df = recalculate_inchikey(salt_df)

    return(salt_df)


def recalculate_inchikey(std_list):
    uncharged_inchikey = []
    for index, row in std_list.iterrows():
        if row['formula_fetched']==row['uncharged_formula']:
            uncharged_inchikey.append(row['inchikey'])
        else:
            uncharged_inchikey.append(Chem.MolToInchiKey(Chem.MolFromSmiles(row['uncharged_smiles'])))
    std_list['uncharged_inchikey']=uncharged_inchikey
    return(std_list)


def complete_smiles(data, reference_col='inchikey', if_na = False):
    if reference_col in ['inchikey', 'name']==False:
        print('allowed values are inchikey or name')
        return

    std_list = data.copy()
    unique_values = data[reference_col].unique()
    smiles_all = []
    for value in tqdm(unique_values):
        smiles_all.append(inchi_to_smiles(value))
    smiles = []
    for index, row in std_list.iterrows():
        idx = list(unique_values).index(row[reference_col])
        smiles.append(smiles_all[idx])
    std_list['smiles_fetched']=smiles
    print('there is '+str(int(std_list['smiles_fetched'].isna().sum()))+' na values in smiles')
    std_list.dropna(subset=['smiles_fetched'], inplace=True)
    return std_list

def complete_formal_charge(data, smiles_column = 'smiles'):
    fcs = []
    for index, row in data.iterrows():
        try:
            mol_temp = Chem.MolFromSmiles(row[smiles_column])
            fcs.append(cal_formal_charge(mol_temp))
        except:
            fcs.append(np.NaN)
    data['formal_charges']=fcs
    return (data)
# def complete_guess_adduct(std_list):
#     data = std_list.copy()
#     guessed_adduct = []
#     mz_offset = []
#     for index, row in std_list_neg.iterrows():
#         guessed_result = guess_adduct(float(row['Precurzor (m/z)']), row['Exact Mass (u)'], '-')
#         guessed_adduct.append(guessed_result[0])
#         mz_offset.append(guessed_result[1])

def guess_adduct(precursor,mono_mass, mode):
    adduct_list = []
    mass_list = []
    for key in single_charged_adduct_mass.keys():
        if key[-1]== mode:
            adduct_list.append(key)
            mass_list.append(single_charged_adduct_mass[key])
    mass_offset = [abs(mono_mass+x-precursor) for x in mass_list]
    idx = np.argmin(mass_offset)
    return(adduct_list[idx], mass_offset[idx])
def complete_mono_mass(std_list, smiles_column = 'smiles', position = 2):
    data = std_list.copy()
    mono_mass = []
    for index, row in data.iterrows():
        try:
            mol_temp = Chem.MolFromSmiles(row[smiles_column])
            mono_mass.append(ExactMolWt(mol_temp))
        except:
            mono_mass.append(np.NaN)
    data.insert(position, "mono_mass", mono_mass)
    return (data)
def complete_formula(data, smiles_column = 'smiles'):
    formulas = []
    for index, row in data.iterrows():
        try:
            mol_temp = Chem.MolFromSmiles(row[smiles_column])
            formulas.append(CalcMolFormula(mol_temp))
        except:
            formulas.append(np.NaN)
    data.insert(2, "Formula_fetched", formulas)
    return (data)
def get_adduct_cols(cols):
    adduct_cols = [x for x in cols if x in single_charged_adduct_mass.keys()]
    return(adduct_cols)


# all_adduct_pos = ['[M]+','[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H2O+H]+']
# all_adduct_neg = ['[M]-','[M-H]-','[M+C2H4O2-H]-','[M-H2O-H]-','[M+FA-H]-','[M+Cl]-','[M+Na-2H]-']
from toolsets.constants import single_charged_adduct_mass
def calculate_precursormz(smiles, adduct):
    mol_ = check_mol(smiles)
    mono_mass =ExactMolWt(mol_)
    if adduct in single_charged_adduct_mass.keys():
        if cal_formal_charge(mol_)==0:
            if adduct not in ['+','-','[M]+','[M+]','[M+]+','[M]-','[M-]','[M-]-']:
                pmz = mono_mass+single_charged_adduct_mass[adduct]
            else:
                pmz = 0
        elif cal_formal_charge(mol_)==1:
            if adduct in ['+','[M]+','[M+]','[M+]+']:
                pmz = mono_mass
            else:
                pmz = 0
        elif cal_formal_charge(mol_)==-1:
            if adduct == ['-','[M]-','[M-]','[M-]-']:
                pmz = mono_mass
            else:
                pmz = 0
    else:
        print(f'the adduct {adduct} you have entered is not a legit adduct')
        return(np.NAN)
    return(np.round(pmz,7))

    # proton = 1.00727646677
    # Na_plus = 22.989218
    # NH4_plus = 18.033823
    # HacH_minus = 59.013851
    # H2OH_minus = -19.01839
    # FaH_minus = 44.998201
    # Cl_minus = 34.969402
    # if cal_formal_charge(mol_)==0:
    #     if(adduct=='[M+NH4]+'):
    #         pmz = ExactMolWt(mol_)+NH4_plus
    #     elif (adduct=='[M+H]+'):
    #         pmz = ExactMolWt(mol_)+proton
    #     elif (adduct=='[M+Na]+'):
    #         pmz = ExactMolWt(mol_)+Na_plus
    #     elif (adduct == '[M+K]+'):
    #         pmz = ExactMolWt(mol_)+single_charged_adduct_mass[adduct]
    #     elif (adduct=='[M-H2O+H]+'):
    #         pmz = ExactMolWt(mol_)-ExactMolWt(Chem.MolFromSmiles('O'))+proton
    #     elif (adduct=='[M-H]-'):
    #         pmz = ExactMolWt(mol_)-proton
    #     elif (adduct=='[M+C2H4O2-H]-'):
    #         pmz = ExactMolWt(mol_)+HacH_minus
    #     elif (adduct=='[M-H2O-H]-'):
    #         pmz = ExactMolWt(mol_)+H2OH_minus
    #     elif (adduct=='[M+FA-H]-'):
    #         pmz = ExactMolWt(mol_)+FaH_minus
    #     elif (adduct=='[M+Cl]-'):
    #         pmz = ExactMolWt(mol_)+Cl_minus
    #     elif (adduct=='[M+Na-2H]-'):
    #         pmz = ExactMolWt(mol_)+20.974666
    #     else:
    #         pmz = 0
    # elif cal_formal_charge(mol_)==1:
    #     if(adduct =='[M]+'):
    #         pmz = ExactMolWt(mol_)
    #     else:
    #         pmz = 0
    # elif cal_formal_charge(mol_)==-1:
    #     if(adduct =='[M]-'):
    #         pmz = ExactMolWt(mol_)
    #     else:
    #         pmz = 0
    # else:
    #     print("you have passed a molecule with multiple charges")
    #     return(np.NAN)


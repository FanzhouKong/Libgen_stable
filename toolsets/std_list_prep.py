import numpy as np
import pandas as pd
from toolsets.helpers import  specify_column
import toolsets.chem_utils as cu
pd.options.mode.chained_assignment = None
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
import os
from toolsets.search import string_search
master_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/Alkaloids_lib/'
filename = 'std_list_raw.csv'
pos_adducts = ['[M+H]+', '[M+Na]+', '[M+NH4]+', '[M]+']
neg_adducts = ['[M-H]-', '[M+Cl]-', '[M+FA-H]-', '[M]-']
identifier = 'smiles'

def std_list_prep(std_list, identifier = 'smiles'):

    col_kept = []

    for keyword in ['name', identifier, 'mix', 'rt_reference']:
        col_name_temp = specify_column(std_list, keyword)
        if len(col_name_temp)>0:
            std_list=std_list.rename(columns={col_name_temp:keyword}, errors="raise")
            col_kept.append(keyword)

    std_list_filtered = std_list[col_kept]

    smiles = []
    mono_mass = []
    formulas = []
    for index, row in std_list_filtered.iterrows():
        desalted_smile = cu.desalter(row['smiles'])
        smiles.append(desalted_smile)
        mono_mass.append(ExactMolWt(Chem.MolFromSmiles(desalted_smile)))
        formulas.append(cu.everything_to_formula(desalted_smile))
    std_list_filtered['smiles']=smiles
    std_list_filtered['formula']=formulas
    std_list_filtered['mono_mass']=mono_mass



    std_list_pos = std_list_filtered.copy()
    std_list_neg = std_list_filtered.copy()

    for adduct in pos_adducts:
        pmzs = []
        for index, row in std_list_pos.iterrows():
            pmzs.append(cu.calculate_precursormz(row['smiles'], adduct, if_smiles=True))
        std_list_pos[adduct]=pmzs
    std_list_pos.mix=std_list_pos.mix.str.replace('neg', 'pos')
    std_list_pos.to_csv(os.path.join(master_dir, 'std_list_pos_curated.csv'), index = False)

    for adduct in neg_adducts:
        pmzs = []
        for index, row in std_list_neg.iterrows():
            pmzs.append(cu.calculate_precursormz(row['smiles'], adduct, if_smiles=True))
        std_list_neg[adduct]=pmzs
    std_list_neg.mix=std_list_neg.mix.str.replace('pos', 'neg')
    return std_list_pos,std_list_neg
    # std_list_neg.to_csv(os.path.join(master_dir, 'std_list_neg_curated.csv'), index = False)
    # std_list_neg.to_csv(os.path.join(master_dir, 'std_list_neg_curated.csv'), index = False)





# below are original, legacy codes
# from toolsets.helpers import find_floats
# import numpy as np
# import toolsets.chem_utils as ch
# def complete_std_list(std_list, smiles_col = 'smiles', adducts = ['[M]+', '[M+H]+', '[M+Na]+','[M+NH4]+']):
#     std_list = complete_formula(std_list, smiles_col = smiles_col)
#     std_list = complete_mono_mass(std_list, smiles_col = smiles_col)
#     std_list = complete_adducts(std_list, smiles_col = smiles_col, adducts = adducts)
#     return(std_list)
# def clean_rt(std_list, colname = 'rt_suggested'):
#     df = std_list.copy()
#     rt_cleaned = []
#     for r in df[colname]:
#         if r == r:
#             rt_cleaned.append(find_floats(r))
#         else:
#             rt_cleaned.append(np.NAN)
#     df[colname]=rt_cleaned
#     return df
# def complete_formula(std_list_raw, smiles_col = 'smiles'):
#     std_list = std_list_raw.copy()
#     formulas = []
#     for index, row in std_list.iterrows():
#         if ch.is_smiles(row[smiles_col])==True:
#             formulas.append(ch.everything_to_formula(row[smiles_col]))
#         else:
#             formulas.append(np.NAN)
#     column_idx = std_list.columns.get_loc(smiles_col)
#     std_list.insert(column_idx+1, 'formula',formulas)
#     return std_list
# def complete_mono_mass(std_list_raw, smiles_col = 'smiles'):
#     std_list = std_list_raw.copy()
#     mono_mass = []
#     for index, row in std_list.iterrows():
#         if ch.is_smiles(row[smiles_col])==True:
#             mono_mass.append(ch.everything_to_mw(row[smiles_col]))
#         else:
#             mono_mass.append(np.NAN)
#     std_list['mono_mass']=mono_mass
#     return std_list
# def complete_adducts(std_list_raw, smiles_col = 'smiles', adducts = '[M+H]+'):
#     std_list = std_list_raw.copy()
#     adducts_dict = {}
#     for adduct in adducts:
#         adducts_dict[adduct]=[]
#     for index, row in std_list.iterrows():
#         if ch.is_smiles(row[smiles_col])==True:
#             for adduct in adducts:
#                 adducts_dict[adduct].append(ch.calculate_precursormz(row[smiles_col], adduct=adduct, if_smiles=True))
#         else:
#             for adduct in adducts:
#                 adducts_dict[adduct].append(np.NAN)
#     for adduct in adducts:
#         std_list[adduct]=adducts_dict[adduct]
#     return std_list
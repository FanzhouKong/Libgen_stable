import os
import sys

# from multiprocess import Pool
import pandas as pd
# from multiprocessing import Pool
from ast import literal_eval
import bisect
import itertools
import time
import toolsets.mass_to_formula as mtf
import numpy as np
from molmass import Formula
from toolsets.search import num_search, string_search, quick_search_sorted
import numpy as np
from toolsets.constants import atom_mass
from molmass import Formula
import toolsets.spectra_operations as so
import toolsets.mass_to_formula as mtf
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import time
from toolsets.std_list_prep import calculate_precursormz
# reference_db_sorted = pd.read_csvl('/Users/fanzhoukong/Documents/GitHub/Libgen_data/formula_db/formulaDB_sorted.csv')
def denoise_dyn(msms, precursor_mz):
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    index_start = np.searchsorted(mass, precursor_mz-1.6,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]
    mass_frag = mass[0:index_start]
    intensity_frag = intensity[0:index_start]
    threshod =np.quantile(intensity_frag, 0.1)
    mass_d = []
    intensity_d = []
    for i in range(0, len(mass_frag)):
        if intensity_frag[i]>threshod:
            mass_d.append(mass_frag[i])
            intensity_d.append(intensity_frag[i])
    if len(mass_d)==0:
        return np.NAN
    mass_d.extend(mass_parent)
    intensity_d.extend(intensity_parent)
    return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)))
def denoise_bp(msms, precursor_mz,threshold = 1):
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    index_start = np.searchsorted(mass, precursor_mz-1.6,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]
    mass_frag = mass[0:index_start]
    intensity_frag = intensity[0:index_start]
    bp_int = np.max(intensity_frag)*threshold/100
    mass_d = []
    intensity_d = []
    for i in range(0, len(mass_frag)):
        if intensity_frag[i]>bp_int:
            mass_d.append(mass_frag[i])
            intensity_d.append(intensity_frag[i])
    if len(mass_d)==0:
        return np.NAN
    mass_d.extend(mass_parent)
    intensity_d.extend(intensity_parent)
    return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)))

    # else:

def denoise_h(msms, smiles, adduct, reference_db_sorted, mass_error = 0.01):
    start = time.time()
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    parent_ion = calculate_precursormz(smiles, adduct)
    # precursor = get_precursor(msms, parent_ion)
    # if precursor == precursor:
    #     parent_ion = precursor
        # state = 'inferred'
        # mass_error = abs(parent_ion-precursor)*3
    # else:
        # state = 'default'
        # mass_error = 0.005
    index_start = np.searchsorted(mass, parent_ion-1.6,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]
    mass_frag = mass[0:index_start]
    intensity_frag = intensity[0:index_start]
    formula = prep_formula(smiles, adduct)
    mass_d = []
    intensity_d = []
    mass_tbd = []
    intensity_tbd = []
    for cur_i in range(len(mass_frag)):
        loss_candidate_white = get_loss_candidate(parent_ion-mass_frag[cur_i], reference_db_sorted, step = mass_error)
        if (check_loss(loss_candidate_white['formstr'].values.tolist(), formula))==True:
            mass_d.append(mass_frag[cur_i])
            intensity_d.append(intensity_frag[cur_i])
        else:
            mass_tbd.append(mass_frag[cur_i])
            intensity_tbd.append(intensity_frag[cur_i])
    check_1 = time.time()
    # all_possible_formulas = get_all_formulas(formula)
    # all_allowed_formulas, all_allowed_masses = get_all_allowed_formula_mass(all_possible_formulas)
    for i in range(len(mass_tbd)):
        loss_mass = parent_ion - mass_tbd[i]
        if check_frag_bl(loss_mass, formula, mass_error):
            mass_d.append(mass_tbd[i])
            intensity_d.append(intensity_tbd[i])
    ei = so.calculate_explained_intensity(so.pack_spectra(mass_d, intensity_d), so.pack_spectra(mass_frag, intensity_frag))
    # print(ei)
    # print(mass_parent)
    if len(mass_d)==0:
        return(np.NAN, 0)
    mass_d.extend(mass_parent)
    intensity_d.extend(intensity_parent)

    # return(mass_d, intensity_d)
    return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)), ei)
def denoise_h_formula(msms, formula, parent_ion, adduct, reference_db_sorted, mass_error = 0.005):
    start = time.time()
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    # parent_ion = calculate_precursormz(smiles, adduct)
    # precursor = get_precursor(msms, parent_ion)
    # if precursor == precursor:
    #     parent_ion = precursor
    # state = 'inferred'
    # mass_error = abs(parent_ion-precursor)*3
    # else:
    # state = 'default'
    # mass_error = 0.005
    index_start = np.searchsorted(mass, parent_ion-1.5,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]
    mass_frag = mass[0:index_start]
    intensity_frag = intensity[0:index_start]
    # formula = prep_formula(smiles, adduct)
    if adduct=='[M+Na]+':
        formula = formula+'Na'
    elif adduct=='[M+NH4]+':
        formula = formula+'NH3'
    elif adduct=='[M+Cl]-':
        formula = formula+'Cl'
    elif adduct=='[M+Hac-H]-':
        formula = formula+'C2H4O2'
    elif adduct =='[M+K]+':
        formula = formula+'K'
    elif adduct == '[M+FA-H]':
        formula = formula +'CH2O2'
    else:
        formula = formula
    formula = (Formula(formula).formula)
    mass_d = []
    intensity_d = []
    mass_tbd = []
    intensity_tbd = []
    for cur_i in range(len(mass_frag)):
        loss_candidate_white = get_loss_candidate(parent_ion-mass_frag[cur_i], reference_db_sorted, step = mass_error)
        if (check_loss(loss_candidate_white['formstr'].values.tolist(), formula))==True:
            mass_d.append(mass_frag[cur_i])
            intensity_d.append(intensity_frag[cur_i])
        else:
            mass_tbd.append(mass_frag[cur_i])
            intensity_tbd.append(intensity_frag[cur_i])
    check_1 = time.time()
    # all_possible_formulas = get_all_formulas(formula)
    # all_allowed_formulas, all_allowed_masses = get_all_allowed_formula_mass(all_possible_formulas)
    for i in range(len(mass_tbd)):
        loss_mass = parent_ion - mass_tbd[i]
        if check_frag_bl(loss_mass, formula, mass_error):
            mass_d.append(mass_tbd[i])
            intensity_d.append(intensity_tbd[i])
    ei = so.calculate_explained_intensity(so.pack_spectra(mass_d, intensity_d), so.pack_spectra(mass_frag, intensity_frag))
    mass_d.extend(mass_parent)
    intensity_d.extend(intensity_parent)
    return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)), ei)
    # return(so.pack_spectra(mass_d, intensity_d), state, mass_error)
def denoise_w(msms, smiles, adduct, reference_db_sorted, mass_error = 0.005):
    start = time.time()
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    parent_ion = calculate_precursormz(smiles, adduct)
    # precursor = get_precursor(msms, parent_ion)
    # if precursor == precursor:
    #     parent_ion = precursor
    # state = 'inferred'
    # mass_error = abs(parent_ion-precursor)*3
    # else:
    # state = 'default'
    # mass_error = 0.005
    index_start = np.searchsorted(mass, parent_ion-1.5,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]
    mass_frag = mass[0:index_start]
    intensity_frag = intensity[0:index_start]
    formula = prep_formula(smiles, adduct)
    mass_d = []
    intensity_d = []
    mass_tbd = []
    intensity_tbd = []
    for cur_i in range(len(mass_frag)):
        loss_candidate_white = get_loss_candidate(parent_ion-mass_frag[cur_i], reference_db_sorted, step = mass_error)
        if (check_loss(loss_candidate_white['formstr'].values.tolist(), formula))==True:
            mass_d.append(mass_frag[cur_i])
            intensity_d.append(intensity_frag[cur_i])
        else:
            mass_tbd.append(mass_frag[cur_i])
            intensity_tbd.append(intensity_frag[cur_i])
    check_1 = time.time()
    # all_possible_formulas = get_all_formulas(formula)
    # all_allowed_formulas, all_allowed_masses = get_all_allowed_formula_mass(all_possible_formulas)
    # for i in range(len(mass_tbd)):
    #     loss_mass = parent_ion - mass_tbd[i]
    #     if check_frag_bl(loss_mass, formula, mass_error):
    #         mass_d.append(mass_tbd[i])
    #         intensity_d.append(intensity_tbd[i])
    ei = so.calculate_explained_intensity(so.pack_spectra(mass_d, intensity_d), so.pack_spectra(mass_frag, intensity_frag))
    mass_d.extend(mass_parent)
    intensity_d.extend(intensity_parent)
    return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)), ei)
def check_frag_bl(frag_mass, formula, mass_error = 0.005):
    nl_candidates = mtf.nl_to_formula(frag_mass, mass_error, formula)
    if nl_candidates != [] and nl_candidates!=['']:
        for nl in nl_candidates:
            if evaluate_nl_blacklist(nl)==True:
                return(True)
        else:
            return(False)
    else:
        return(False)
# if intensity_precursor[0] != -1:
#     mass_d.extend(mass_precursor)
#     intensity_d.extend(intensity_precursor)

def find_loss(loss_mass, mass_error, all_allowed_formulas, all_allowed_masses):
    index_start = np.searchsorted(all_allowed_masses, loss_mass-mass_error,side = 'left')
    index_end = np.searchsorted(all_allowed_masses, loss_mass+mass_error,side = 'right')
    if len(all_allowed_masses[index_start:index_end])>0:
        return(all_allowed_formulas[index_start:index_end])
    else:
        return ([])
def get_all_allowed_formula_mass(all_possible_formulas):
    all_allowed_formulas = []
    all_allowed_masses = []
    check_21 = time.time()
    for formula in all_possible_formulas:
        if evaluate_nl_blacklist(formula)==True:
            all_allowed_formulas.append(formula)
            all_allowed_masses.append(Formula(formula).isotope.mass)
    if len(all_allowed_formulas)<1:
        return([], [])
    allowed_mass_sorted, allowed_formula_sorted = zip(*sorted(zip(all_allowed_masses, all_allowed_formulas)))
    return(allowed_formula_sorted,allowed_mass_sorted )
def get_all_formulas(formula):
    precursor_formula = MolecularFormula()
    precursor_formula.from_string(Formula(formula).formula)
    precursor_data = precursor_formula.get_data()
    formula_range = [range(x + 1) for x in precursor_data]
    all_possible_candidate_formula = np.array(
        list(itertools.product(*formula_range)), numpy_formula_format)
    all_possible_formulas = translate_matrix_into_formula(all_possible_candidate_formula)
    return(all_possible_formulas)
def translate_matrix_into_formula(matrix):
    formulas = []
    for i in matrix:
        if np.max(i)>0:
            formula_temp = ''
            for j in range(len(i)):
                if i[j]>0:
                    if i[j]>1:
                        formula_temp = formula_temp+atom_list[j]+str(i[j])
                    elif i[j]==1:
                        formula_temp = formula_temp+atom_list[j]
            formulas.append(formula_temp)
    return formulas
def get_precursor(msms, precursor_mz_tho):
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    index_start = np.searchsorted(mass, precursor_mz_tho-1.5,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]
    if len(mass_parent)>0:
        precursor_idx = np.argmax(intensity_parent)
        return(mass_parent[precursor_idx])
    else:
        return(np.NAN)

from rdkit.Chem.AllChem import CalcExactMolWt
from rdkit import Chem
from toolsets.constants import single_charged_adduct_mass
def denoise_hybrid_msp(instance, reference_db_sorted, typeofmsms = 'peaks_recalibrated', mass_error =0.02 ):
    mass, intensity = so.break_spectra(instance[typeofmsms])
    adduct = instance['Precursor_type']
    # parention = float(instance['PrecursorMZ'])
    parention = CalcExactMolWt(Chem.MolFromSmiles(instance['SMILES']))+single_charged_adduct_mass[adduct]
    
    formula = prep_formula(instance['SMILES'], adduct)
    mass_d = []
    intensity_d = []
    for cur_i in range(len(mass)):
        loss_candidate_white = get_loss_candidate(parention-mass[cur_i], reference_db_sorted, step = mass_error)
        if (check_loss(loss_candidate_white['formstr'].values.tolist(), formula))==True:
            mass_d.append(mass[cur_i])
            intensity_d.append(intensity[cur_i])
        else:
            loss_candidate_black = mtf.nl_to_formula(parention-mass[cur_i], mass_error, formula)
            if(check_candidate_black(loss_candidate_black))==True:
                mass_d.append(mass[cur_i])
                intensity_d.append(intensity[cur_i])
    return(so.pack_spectra(mass_d, intensity_d))



def check_candidate_black(loss_candidate_black):
    if len(loss_candidate_black) > 0 and loss_candidate_black!=['']:
        for nl in loss_candidate_black:
            if evaluate_nl_blacklist(nl)==True:
                return(True)
        return(False)
    else:
        return(False)
def get_loss_candidate(mass, reference_db, step):
    mass_neg = mass- atom_mass['H']
    mass_pos = mass+atom_mass['H']
    losses = pd.DataFrame()
    for m in [mass_neg, mass, mass_pos]:
        # print(m)
        losses = losses.append(quick_search_sorted(reference_db, 'mass', m, step))
    return (losses)
def check_loss(loss_formula_candidate, formula):
    # print(formula)
    formula_dict = parse_formula_dict(Formula(formula))
    for loss_formula in loss_formula_candidate:
        loss_temp_dict = parse_formula_dict(Formula(loss_formula))
        if check_candidate(loss_temp_dict, formula_dict)==True:
            # print(loss_formula)
            return(True)
    return(False)
def check_candidate(loss_temp_dict, formula_dict):
#     both loss_temp_d and formula_d are in dataframe form
    if set(loss_temp_dict.keys()).issubset(set(formula_dict.keys())):
        # is all loss atom subset of formula atom?
        for element in loss_temp_dict.keys():
            if loss_temp_dict[element] > formula_dict[element]:
                return( False)
        return(True)
    else:
        return(False)
def denoise_blacklist(instance, typeofmsms = 'peaks_recalibrated', mass_error =0.02):
    # print("i am in new")
    import toolsets.spectra_operations as so
    mass, intensity = so.break_spectra(instance[typeofmsms])
    parention = float(instance['reference_precursor_mz'])
    # pep_mass, pep_intensity = get_parentpeak(mass, intensity, parention)
    adduct = instance['reference_adduct']
    formula = prep_formula(instance['reference_smiles'], adduct)
    # mass_frag, intensity_frag, mass_precursor, intensity_precursor = remove_precursor(mass, intensity, parention)

    # if instance['reference_adduct']=='[M-H2O-H]-' or instance['reference_adduct']=='[M+H2O+H]+':
    #     formula = formula+'H2O'
    mass_d = []
    intensity_d = []

    # print(threshold)
    for i in range(0,len(mass)):
        # print("i am evaluating", mass[i])
        nl_candidates = mtf.nl_to_formula(parention - mass[i], mass_error, formula)
        if nl_candidates != [] and nl_candidates!=['']:
            for nl in nl_candidates:
                if evaluate_nl_blacklist(nl)==True:
                    mass_d.append(mass[i])
                    intensity_d.append(intensity[i])
                    break
    # if intensity_precursor[0] != -1:
    #     mass_d.extend(mass_precursor)
    #     intensity_d.extend(intensity_precursor)
    return(so.pack_spectra(mass_d, intensity_d))
def denoise_single(frag_ion, parention, smiles, mass_error = 0.02, adduct = '[M+H]+',ifppm = False, ifnl = True):
    import toolsets.spectra_operations as so
    parention = float(parention)
    threshold = so.set_tolerance(mass_error = mass_error, ifppm=ifppm, precursormz=parention)
    # formula= str(formula)
    formula = prep_formula(smiles, adduct)
    if ifnl ==True:
        nl_candidates = mtf.nl_to_formula(parention - frag_ion, threshold, formula)
        return(nl_candidates)
    if ifnl == False:
        nl_candidates = mtf.nl_to_formula(frag_ion, threshold, formula)
        return(nl_candidates)

def prep_formula(smiles, adduct):
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    benzene_pattern = Chem.MolFromSmiles('C1=CC=CC=C1')
    if mol.HasSubstructMatch(benzene_pattern) ==True:
        formula = CalcMolFormula(mol)
        formula = formula +'N2'
    else:
        formula = CalcMolFormula(mol)
    # formula = CalcMolFormula(mol)

    # formula= str(formula)
    if adduct=='[M+Na]+':
        formula = formula+'Na'
    elif adduct=='[M+NH4]+':
        formula = formula+'NH3'
    elif adduct=='[M+Cl]-':
        formula = formula+'Cl'
    elif adduct=='[M+Hac-H]-':
        formula = formula+'C2H4O2'
    elif adduct =='[M+K]+':
        formula = formula+'K'
    elif adduct == '[M+FA-H]':
        formula = formula +'CH2O2'
    elif formula[-1]=='+' or formula[-1]=='-':
        formula = formula[:-1]
    else:
        formula = formula
    return(Formula(formula).formula)
def post_processing(data, raw_column = 'peaks_recalibrated',denoised_column = 'peaks_denoised'):
    import toolsets.spectra_operations as so
    ei = []
    for index, row in data.iterrows():
        ei.append(so.calculate_explained_intensity(row[raw_column], row[denoised_column]))
    data['ei']=ei
    peaks_denoisied_normalized = []
    normalized_entropy = []
    spectrum_entropy = []
    # data['c_id'] = np.arange(len(data))
    for index, row in data.iterrows():
        peaks_denoisied_normalized.append(so.normalize_spectrum(row[denoised_column]))
        normalized_entropy.append(so.normalized_entropy(row[denoised_column], order=4))
        spectrum_entropy.append(so.spectral_entropy(row[denoised_column]))
    data['peaks_denoised_normalized']=peaks_denoisied_normalized
    data['spectrum_entropy']=spectrum_entropy
    data['normalized_entropy']=normalized_entropy
    return(data)


from toolsets.mass_to_formula import atom_list
from toolsets.mass_to_formula import MolecularFormula, atom_mass_array, numpy_formula_format,atom_list
def get_allowed_formula(all_possible_candidate_formula, element):

    idx = atom_list.index(element)
    all_allowed_candidate_formula = np.copy(all_possible_candidate_formula)
    pop_idx = []
    for i in range(len(all_possible_candidate_formula)):
        if np.max(all_possible_candidate_formula[i][0:idx])==0 and np.max(all_possible_candidate_formula[i][idx+1:])==0:
            pop_idx.append(i)
    all_allowed_candidate_formula = np.delete(all_allowed_candidate_formula, pop_idx, axis=0)
    # all_allowed_candidate_formula.append(all_possible_candidate_formula[i])
    return(all_allowed_candidate_formula)
def get_all_allowed_formula(all_possible_candidate_formula, elements):
    all_allowed_candidate_formula = np.copy(all_possible_candidate_formula)
    for element in elements:
        all_allowed_candidate_formula = get_allowed_formula(all_allowed_candidate_formula, element = element)
    return(all_allowed_candidate_formula)
def math_check(smiles, adduct, return_mass_list = False):
    precursor_mz = calculate_precursormz(smiles, adduct)
    formula = prep_formula(smiles, adduct)
    precursor_formula = MolecularFormula()

    precursor_formula.from_string(Formula(formula).formula)
    precursor_data = precursor_formula.get_data()
    formula_range = [range(x + 1) for x in precursor_data]
    forbidden_elements = ['C','N']
    if 'Na' in formula:
        forbidden_elements.append('Na')
    elif 'K' in formula:
        forbidden_elements.append('K')
    elif 'Cl' in formula:
        forbidden_elements.append('Cl')
    all_possible_candidate_formula = np.array(
        list(itertools.product(*formula_range)), numpy_formula_format)
    all_allowed_candidate_formula = get_all_allowed_formula(all_possible_candidate_formula, elements=forbidden_elements)
    all_possible_mass = np.sum(
        atom_mass_array * all_allowed_candidate_formula, axis=1)
    all_possible_mass.sort()
    diff = np.diff(all_possible_mass)
    allowed_length = len(diff[diff>0.005])*0.005*2
    # if return_mass_list == True:
    #     return(all_possible_mass)
    # peak_list = make_peak_list(all_possible_mass)
    # peak_list_checked = check(peak_list)
    # allowed_region = 0
    # for i in range(0, len(peak_list_checked)):
    #     allowed_region = allowed_region+(peak_list_checked[i][2]-peak_list_checked[i][0])
    return(allowed_length/precursor_mz*100)
def make_peak_list(all_mass_list):
    peak_list = []
    all_mass_list.sort()
    for mass in all_mass_list:
        peak_list.append([mass-0.005, mass, mass+0.005])
    return peak_list
def check_overlap(peak_list):
    # peak_list.sort()
    for p in range(0, len(peak_list)-1):
        if peak_list[p][2]>peak_list[p+1][0]:
            return(p)
    return(-1)
def check(peak_list):

    while check_overlap(peak_list)!= -1:
        idx = check_overlap(peak_list)
        peak_list[idx][2]=peak_list[idx+1][2]
        peak_list = np.delete(peak_list, idx+1, axis=0)

    return(peak_list)
    # peak_list[p][2]=peak_list[p+1][2]
    # data_good=num_search(data, 'ei', ei_threshold, direction=">", inclusion=True)
    #     # data=num_search(data, 'ei', ei_threshold, direction="<", inclusion=False)
    # data_good =num_search(data_good, 'spectrum_entropy', 0.5, direction='>', inclusion=True)
    # data_good.reset_index(drop = True, inplace= True)
    # if high_quality ==True:
    #     return(data_good)
    # else:
    #     # data['c_id'] = np.arange(len(data_HCD))
    #     data_bad = data[~data['library_id'].isin(data_good['library_id'])]
    #     # data_bad =num_search(data_bad, 'spectrum_entropy', 0.5, direction='>', inclusion=True)
    #     data_bad.reset_index(drop = True, inplace= True)
    #
    #     return(data_bad)
# def prep_formula(formula, adduct = ['[M+H]+']):
#     formula= str(formula)
#     if adduct=='[M+Na]+':
#         formula = formula+'Na'
#     if adduct=='[M+NH4+]':
#         formula = formula+'NH3'
#     if adduct=='[M+Cl]-':
#         formula = formula+'Cl'
#     if adduct=='[M+C2H4O2-H]-':
#         formula = formula+'C2H4O2'
#     if formula[-1]=='+' or formula[-1]=='-':
#         formula = formula[:-1]
#     return(formula)
# def find_parention(mass, intensity, precursormz, mass_error = 10, ifppm = True):
#     precursormz = float(precursormz)
#     tol = so.set_tolerance(mass_error, ifppm, precursormz)
#     # mass_raw, intensity_raw = so.break_spectra(msms)
#     # if precursormz <=200:
#     #     tolerance = 0.01
#     # if precursormz >200:
#     #     tolerance = precursormz*tolerance/1E6
#
#     # mass_sorted, intensity_sorted = zip(*sorted(zip(mass_raw, intensity_raw)))
#     lower_bound = precursormz-tol
#     upper_bound = precursormz+tol
#     # print(" i am here!")
#     lower_bound_i = bisect.bisect_left(mass, lower_bound)
#     # print("i have passed lower")
#     upper_bound_i = bisect.bisect_right(mass, upper_bound, lo=lower_bound_i)
#     mass_candi = mass[lower_bound_i:upper_bound_i]
#     intensity_candi = intensity[lower_bound_i:upper_bound_i]
#     # return(mass_candi)
#     if mass_candi ==[]:
#         return(precursormz)
#     else:
#         max_index = intensity_candi.index(max(intensity_candi))
#         return(mass_candi[max_index])

    # return (mass_candi, max_index)
# def get_parentpeak(mass, intensity, parention, mass_error = 0.01, ifppm = False):
#     parention = float(parention)
#     tolerance = so.set_tolerance(mass_error, ifppm, parention)
#     lower_bound = parention-tolerance
#     upper_bound = parention+tolerance

#     lower_bound_i = bisect.bisect_left(mass, lower_bound)
#     upper_bound_i = bisect.bisect_right(mass, upper_bound, lo=lower_bound_i)

#     mass_candi = mass[lower_bound_i:upper_bound_i]
#     # print("checkpoint!")
#     intensity_candi = intensity[lower_bound_i:upper_bound_i]
#     # return(mass_candi)
#     if mass_candi !=[]:
#         max_index = intensity_candi.index(max(intensity_candi))
#         return(mass_candi[max_index], intensity_candi[max_index])
#     else:
#         return(parention, 0)
def parse_formula_dict(f):
    f = f.composition()
    element =[]
    count = []
    for i in range(len(f)):
        element.append(f[i][0])
        count.append(f[i][1])
    formula_dict = dict(zip(element, count))
    return(formula_dict)
def evaluate_nl_blacklist(nl):
    if len(nl)==0:
        return False

    mol = MolecularFormula()
    mol.from_string(nl)
    if np.max(np.append(mol.get_data()[0:1], mol.get_data()[2:]))==0 or np.max(np.append(mol.get_data()[0:2], mol.get_data()[3:]))==0:

        return False
    n_H = mol.get_data()[atom_list.index('H')]
    n_C = mol.get_data()[atom_list.index('C')]
    n_N = mol.get_data()[atom_list.index('N')]
    n_O = mol.get_data()[atom_list.index('O')]
    n_F = mol.get_data()[atom_list.index('F')]
    n_Na = mol.get_data()[atom_list.index('Na')]
    n_P = mol.get_data()[atom_list.index('P')]
    n_S = mol.get_data()[atom_list.index('S')]
    n_Cl = mol.get_data()[atom_list.index('Cl')]
    n_K = mol.get_data()[atom_list.index('K')]
    n_Br = mol.get_data()[atom_list.index('Br')]
    n_I = mol.get_data()[atom_list.index('I')]
    if n_C != 0 or n_O != 0 or n_N != 0 or n_P != 0 or n_S != 0:
        if n_H+n_Br+n_I+n_F+n_Cl>n_C*4+n_O*4+n_N*5+n_P*5+n_S*4:

            return False
    if n_Na >1:
        return(False)
    if n_N/n_C >4:
        return(False)
    if n_O/n_C>3:
        return(False)
    if n_H/n_C>9:
        return(False)
    if n_P/n_C>2:
        return False
    if n_K >1:
        return(False)
    if n_Na >0 and n_H >0:
        if np.max(np.append(mol.get_data()[1:5], mol.get_data()[6:]))==0:
            return False
    if n_K >0 and n_H >0:
        if np.max(np.append(mol.get_data()[1:9], mol.get_data()[10:]))==0:
            return False
    return True

# def evaluate_nl_blacklist(nl):
#     if len(nl)==0:
#         return(False)
#     if len(Formula(nl).composition())==1:
#         if Formula(nl).composition()[0][0] =='C' or Formula(nl).composition()[0][0] =='N':
#             return(False)
#
#         else:
#             return(True)
#
#     else:
#         #
#         # if:
#         #     pass
#         # else:
#         return(True)



# def evaluate_frag(instance):
#     allowed_formula = []
#     precursormz = float(spec['PrecursorMZ'])
#     mass, intensity = prep_msms(spec['spectrum'], ifnist=True)
#     parention = find_parention(mass, intensity, precursormz)
#     mass, intensity = truncate_msms(mass, intensity, parention)
#     for i in range(0, len(mass)):
#         allowed_formula_temp = mtf.nl_to_formula(parention - mass[i], 10, spec['Formula'])
#         if allowed_formula_temp!=[] and allowed_formula_temp!=['']:
#             allowed_formula.extend(allowed_formula_temp)
#             # return(allowed_formula)
#         # else:
#         #     allowed_formula.extend([-1])
#             # return(allowed_formula)
#     return (allowed_formula)
def find_losses_nist(spec):
    allowed_formula = []
    precursormz = float(spec['PrecursorMZ'])
    mass, intensity = prep_msms(spec['spectrum'], ifnist=True)
    parention = get_parentpeak(mass, intensity, precursormz)
    mass, intensity = truncate_msms(mass, intensity, parention)
    for i in range(0, len(mass)):
        allowed_formula_temp = mtf.nl_to_formula(parention - mass[i], 10, spec['Formula'])
        if allowed_formula_temp!=[] and allowed_formula_temp!=['']:
            allowed_formula.extend(allowed_formula_temp)
            # return(allowed_formula)
        # else:
        #     allowed_formula.extend([-1])
            # return(allowed_formula)
    return (allowed_formula)





    # pep_mass, pep_intensity = get_parentpeak(mass, intensity, precursormz)
    # mass, intensity = truncate_msms(mass, intensity, pep_mass)
    # if adduct=='[M+Na]+':
    #     formula = formula+'Na'
    # if adduct=='[M+NH4+]':
    #     formula = formula+'NH3'
    # if adduct=='[M+Cl]-':
    #     formula = formula+'Cl'
    # mass_d = []
    # intensity_d = []
    # threshold = so.set_tolerance(mass_error, ifppm = ifppm, precursormz=precursormz)
    #
    # # print(threshold)
    # for i in range(0,len(mass)):
    #     # print("i am evaluating", mass[i])
    #     nl_candidates = mtf.nl_to_formula(pep_mass - mass[i], threshold, formula)
    #     if nl_candidates != [] and nl_candidates!=['']:
    #         for nl in nl_candidates:
    #             if evaluate_nl_blacklist(nl)==True:
    #                 mass_d.append(mass[i])
    #                 intensity_d.append(intensity[i])
    #                 break
    # if pep_intensity != 0:
    #     mass_d.append(pep_mass)
    #     intensity_d.append(pep_intensity)
    # return(so.sort_spectra(so.pack_spectra(mass_d, intensity_d)) )
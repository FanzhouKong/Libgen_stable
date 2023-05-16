#!/bin/env/python3
import itertools
import numpy as np
import pandas as pd
import re
import math
import sys
from molmass import Formula
from toolsets.constants import single_charged_adduct_mass, atom_valence, atom_mass
import toolsets.denoising_related_functions as de
# from toolsets.denoising_related_functions import evaluate_nl_blacklist
from brainpy import isotopic_variants
from toolsets.search import string_search
import toolsets.spectra_operations as so
# from toolsets.spectra_operations import normalize_spectrum, break_spectra, pack_spectra, set_tolerance, standardize_spectra, num_peaks
from sklearn.linear_model import LinearRegression
atom_list = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']
numpy_formula_format = np.int16
# M, f -> measured, m, p-> therotical
def calculate_isotope_simple(ms1_silico, ms1):
    ms1_silico =standardize_spectra(ms1_silico)
    ms1 = standardize_spectra(ms1)
    mass_silico, intensity_silico = break_spectra(ms1_silico)
    mass, intensity = break_spectra(ms1)
    score= 1
    for cur_i in range(0, num_peaks(ms1_silico)):
        score = score*(1-abs(intensity[cur_i]-intensity_silico[cur_i]))
    return score

def get_spec_p_intensity(ms1_silico, ms1, off_set = 0.002):
    mj, pj = break_spectra(ms1_silico)
    Mj,  fj = break_spectra(ms1)
    fj = [item +off_set for item in fj]
    fj = [item/sum(fj) for item in fj]
    p_intensity = 1
    for cur_i in range(0, len(fj)):
        # print(get_p_intensity(fj[cur_i], pj[cur_i]))
        p_intensity = p_intensity*get_p_intensity(fj[cur_i], pj[cur_i])
        # print(p_intensity)
    return(p_intensity)
def get_spec_p_mass(ms1_silico, ms1, include_base = False):
    mj, pj = break_spectra(ms1_silico)
    Mj,  fj = break_spectra(ms1)
    M0 = Mj[0]
    m0 = mj[0]
    p_mass = 1
    if include_base ==True:
        start = 0
    else:
        start = 1
    for cur_i in range(start, len(Mj)):
        if cur_i ==0:
            p_mass = p_mass*get_p_mass(Mj[cur_i], mj[cur_i], fj[cur_i])
        else:      
            p_mass = p_mass*get_p_mass(Mj[cur_i], mj[cur_i], fj[cur_i], M0 = M0, m0 = m0)         
    return(p_mass)
def get_p_intensity(fj, pj):
    # fj is the intensity of measured peak
    sigma = get_sigma_intensity(fj)
    p = math.erfc(abs(math.log( fj/pj, 10))/(math.sqrt(2)*sigma))
    return(p)

def get_p_mass(Mj, mj, fj, M0 = None, m0 = None):
    sigma = get_sigma_mass(fj)
    # print(sigma)
    # print('alpha1, alpha 0',alpha1, alpha0)
    if M0 !=None and m0!=None:
        # calculate by difference to eliminate sysmtematic error
        Mj = Mj-M0
        mj = mj-m0
    #     this is when comaring monoisotopic peaks
    p = math.erfc(abs(Mj-mj)/(math.sqrt(2)*sigma))
    return(p)
def get_sigma_mass(intensity, mass_accuray = 0.005):
    alpha = get_alpha(intensity)
    return(1/3*alpha*mass_accuray)
def get_sigma_intensity(intensity):
    beta = get_beta(intensity)
    return(1/3*math.log((1+beta),10 ))

def get_beta(intensity):
    if intensity >1:
        intensity = intensity/100
    if intensity >=0.8:
        return (0.08)
    elif intensity >=0.5 and intensity < 0.8:
        x = np.array([0.5, 0.8]).reshape((-1, 1))
        y = np.array([0.12, 0.08])

    elif intensity >=0.2 and intensity <0.5:
        x = np.array([0.2, 0.5]).reshape((-1, 1))
        y = np.array([0.25, 0.12])
    elif intensity >=0.1 and intensity <0.2:
        x = np.array([0.1,0.2 ]).reshape((-1, 1))
        y = np.array([0.5, 0.25])
    elif intensity >=0.05 and intensity <0.1:
        x = np.array([0.05, 0.1]).reshape((-1, 1))
        y = np.array([1.2, 0.5])
    elif intensity >=0.01 and intensity <0.05:
        x = np.array([0.01, 0.05]).reshape((-1, 1))
        y = np.array([2, 1.2])
    elif intensity <0.01:
        return(2)
    model = LinearRegression().fit(x, y)
    return (model.predict(np.array([intensity]).reshape(1,-1))[0])
def get_alpha(intensity):
    if intensity >1:
        intensity = intensity/100
    if intensity >=0.8:
        return (1)
    elif intensity >=0.5 and intensity < 0.8:
        x = np.array([0.8, 0.5]).reshape((-1, 1))
        y = np.array([1, 1.1])
    elif intensity >=0.2 and intensity <0.5:
        x = np.array([0.5, 0.2]).reshape((-1, 1))
        y = np.array([1.1, 1.2])
    elif intensity >=0.1 and intensity <0.2:
        x = np.array([0.2, 0.1]).reshape((-1, 1))
        y = np.array([1.2, 2])
    elif intensity >=0.05 and intensity <0.1:
        x = np.array([0.05, 0.1]).reshape((-1, 1))
        y = np.array([2.5, 2])
    elif intensity >=0.01 and intensity <0.05:
        x = np.array([0.01, 0.05]).reshape((-1, 1))
        y = np.array([2.5, 10])
    elif intensity <0.01:
        return(10)
    model = LinearRegression().fit(x, y)
    return (model.predict(np.array([intensity]).reshape(1,-1))[0])

def get_candidate(mass, mass_error, ifppm):
    result = neutral_to_formula(mass, 
        set_tolerance(mass_error=mass_error,  ifppm=ifppm, precursormz=mass))
    compliance = []
    formulas = []
    for i in (range(len(result))):
        compliance.append(check_chemistry(result[i].__str__()))
        # f_temp = Formula
        formulas.append(Formula(result[i].__str__()).formula)
        # RDBE.append(get_RDBE(Formula(result[i].__str__())))
    tentative_result = pd.DataFrame(list(zip(formulas, compliance )), columns=['Formulas', 'chemistry'])
    tentative_result = string_search(tentative_result, 'chemistry', True)
    tentative_result.reset_index(inplace = True, drop = True)
    return(tentative_result)
atom_dict = {a: i for i, a in enumerate(atom_list)}
len_atom_dict = len(atom_dict)
atom_mass_array = np.zeros(len_atom_dict, np.float32)
for atom in atom_dict:
    atom_mass_array[atom_dict[atom]] = atom_mass[atom]


def neutral_to_formula(mass, mass_error):
    mass = float(mass)
    mass_error = float(mass_error)
    return precursor_mass_to_formula(mass+1.007276, mass_error, '[M+H]+')
def parse_formula_dict(f):
    f = f.composition()
    element =[]
    count = []
    for i in range(len(f)):
        element.append(f[i][0])
        count.append(f[i][1])
    formula_dict = dict(zip(element, count))
    return(formula_dict)
def parse_formula(f):
    f = f.composition()
    element =[]
    count = []
    for i in range(len(f)):
        element.append(f[i][0])
        count.append(f[i][1])
    df = pd.DataFrame(list(zip(element, count)),
               columns =['element', 'count'])
    return(df)
def get_nvalence(formula):
    # this formula is already a molmass.Formula item!!!!
    df = parse_formula(formula)
    total_valence = 0
    # return(df)
    for index, row in df.iterrows():
        total_valence=total_valence+atom_valence[row['element']]*row['count']
    return(total_valence)

def get_num_atom(df, atoms):

    # atoms has to be in list
    df_atoms = df[df['element'].isin(atoms)]
    return(df_atoms['count'].sum())

def get_RDBE(formula):
    # formula = Formula(formula)
    df = parse_formula(formula)
    RDBE = get_num_atom(df, ['C'])-(get_num_atom(df, ['H'])+get_num_atom(df, ['F','Cl','Br','I']))/2+get_num_atom(df, ['N'])/2+1
    return RDBE

def check_senior3(formula):
    # print("i am in new")

    # f = Formula(formula)
    # RDBE = get_RDBE(f)
    natoms = formula.atoms
    nvalence = get_nvalence(formula)
    if nvalence >= 2*natoms-1:
        return (True)
    else:
        return (False)

def check_chemistry(formula):
    formula = Formula(formula)
    nvalence = get_nvalence(formula)
    if check_senior3(formula) == True and get_RDBE(formula) >-0.5 and de.evaluate_nl_blacklist(formula.formula) and nvalence%2==0:
        return (True)
    else:
        return (False)

from brainpy import isotopic_variants
# from toolsets.constants import single_charged_adduct_mass
def get_isotope(formula, digit = 3, adduct = None):
    formula_temp = Formula(formula)
    formula_temp = parse_formula(formula_temp)
    mol_temp = dict(zip(formula_temp['element'], formula_temp['count']))
    theoretical_isotopic_cluster = isotopic_variants(mol_temp, npeaks=digit, charge=0)
    mass = []
    intensity = []
    for i in range(digit):
        if adduct != None:
            mass.append(theoretical_isotopic_cluster[i].mz+single_charged_adduct_mass[adduct])
        else:
            mass.append(theoretical_isotopic_cluster[i].mz)
        intensity.append(theoretical_isotopic_cluster[i].intensity)
    ms1 = pack_spectra(mass, intensity)
    ms1 = normalize_spectrum(ms1)
    return (ms1)

def _calculate_formula(mass_start, mass_end, candidate_formula_array, cur_i, result):
    # print(cur_i)
    atom_mass_cur = atom_mass_array[cur_i]
    # print(atom_mass_cur)
    atom_num = math.floor(mass_end / atom_mass_cur)
    if cur_i == 0:
        # print('i am in if loop')
        # This is H
        h_num_low = mass_start / atom_mass_cur
        if atom_num >= h_num_low:
            # print("i am in ifif")
            candidate_formula_array[0] = atom_num
            result.append(MolecularFormula(candidate_formula_array))
        # else:
        #     candidate_formula_array[0] = atom_num-1
        #     result.append(MolecularFormula(candidate_formula_array))
    else:
        # print('i am in else loop')
        if atom_num>0:
            for i in range(atom_num+1):
                f = np.copy(candidate_formula_array)
                f[cur_i] = i
                # print(mass_start - i * atom_mass_cur)
                # print(atom_list[cur_i])
                _calculate_formula(mass_start - i * atom_mass_cur, mass_end - i * atom_mass_cur,
                                   f, cur_i - 1, result)

        else:
            f = np.copy(candidate_formula_array)
            f[cur_i] = 0
            _calculate_formula(mass_start, mass_end,
                                   f, cur_i - 1, result)

def precursor_mass_to_formula(mass, mass_error, addition):
    mol_mass = mass - single_charged_adduct_mass[addition]
    # print(mol_mass)
    lo_mass = mol_mass - mass_error
    hi_mass = mol_mass + mass_error
    # print(mass_error)

    result = []
    candidate_formula_array = np.zeros(len_atom_dict, numpy_formula_format)
    _calculate_formula(lo_mass, hi_mass, candidate_formula_array,
                       len(candidate_formula_array) - 1, result)
    return result






class MolecularFormula(object):
    __slots__ = ['_data', '_hash']

    def __init__(self, data=None):
        if data is not None:
            self._data = np.array(
                data, numpy_formula_format, copy=True)
        else:
            self._data = np.zeros(len_atom_dict, numpy_formula_format)
        self._hash = None
        pass

    def __getitem__(self, item):
        return self._data[atom_dict[item]]

    def __setitem__(self, key, value):
        self._data[atom_dict[key]] = value

    def __str__(self):
        string = ''

        for atom in atom_list:
            atom_num = self[atom]
            if atom_num:
                if atom_num > 1:
                    string += atom + str(atom_num)
                else:
                    string += atom
        return string

    def from_string(self, key_string):
        f = Formula(key_string)
        # f.composition()
        all_atom_nums = []
        for i in range(len(f.composition())):
            all_atom_nums.append((f.composition()[i][0], f.composition()[i][1]))
        for atom_num in all_atom_nums:
            self[atom_num[0]] = int(atom_num[1])
        pass

    def get_data(self):
        return self._data

    def get_mass(self):
        return np.sum(atom_mass_array * self._data)

# def product_mass_to_formula(mass, precursor_formula, mass_error = 0.01, adduct=None ):
#     if adduct != None:
#         mass -= single_charged_adduct_mass[adduct]

#     lo_mass = mass - mass_error
#     hi_mass = mass + mass_error

#     # Generate candidate range
#     precursor_data = MolecularFormula()
#     precursor_data.from_string(precursor_formula)
#     formula_range = [range(x + 1) for x in precursor_data.get_data()]
#     all_possible_candidate_formula = np.array(
#         list(itertools.product(*formula_range)), numpy_formula_format)
#     all_possible_mass = np.sum(
#         atom_mass_array * all_possible_candidate_formula, axis=1)

#     candidate_data = all_possible_candidate_formula[(lo_mass <= all_possible_mass) & (all_possible_mass <= hi_mass)]

#     result = []
#     for data in candidate_data:
#         formula = MolecularFormula(data)
#         result.append(formula)
#     return result
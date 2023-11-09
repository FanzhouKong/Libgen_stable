import requests
import numpy as np
import pandas as pd
from tqdm import  tqdm
from urllib.request import urlopen
from urllib.parse import quote

def name_to_smiles(ids):
    try:
        r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ids}/property/CanonicalSMILES/JSON').json()
        return r['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except:
        return np.NAN
def smiles_to_name(smile):
    try:
        r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smile}/property/IUPACName/JSON').json()
        return r['PropertyTable']['Properties'][0]['IUPACName']
    except:
        return(np.NAN)
def inchi_to_smiles(inchikey=None , show_example = False):
    inputt='InChIKey'
    outputt ='CanonicalSMILES'
    full_inchi = inchikey
    half_inchi = inchikey[0:14]
    if show_example == True:
        print("allowed inputt:InChIKey, CanonicalSMILES, IsomericSMILES, MolecularFormula, XLogP...")
        return

    r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{inputt}/{full_inchi}/property/{outputt}/JSON').json()
    if len(r.keys())==1 and list(r.keys())[0] == 'Fault': #check if full length is not found
        smiles_gnps=GNPS_get(inputt = inputt.lower(), content = full_inchi, outputt = outputt[-6:].lower(),firstblock_only = False)
        # return(smiles_gnps)
        if smiles_gnps == smiles_gnps:
            return(smiles_gnps)
        else:
            r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{inputt}/{half_inchi}/property/{outputt}/JSON').json()
            # return(r)
            if len(r.keys())==1 and list(r.keys())[0] == 'Fault':
                smiles_gnps=GNPS_get(inputt = inputt.lower(), content = full_inchi, outputt = outputt[-6:].lower(),firstblock_only = True)
                if smiles_gnps == smiles_gnps:
                    return smiles_gnps
                else:
                    return(np.NAN)
            else:
                neutral_one = get_neutral_one(r)
                return(neutral_one)
    else:# smiles found by full inchikey!
        neutral_one = get_neutral_one(r)
        return(neutral_one)
def get_neutral_one(r):
    from toolsets.std_list_prep import cal_formal_charge
    if len(r['PropertyTable']['Properties'])>1:
        abs_formal_charge = []
        for t in r['PropertyTable']['Properties']:
            abs_formal_charge.append(abs(cal_formal_charge(t['CanonicalSMILES'])))
        return_idx = np.argmin(abs_formal_charge)
        return(r['PropertyTable']['Properties'][return_idx]['CanonicalSMILES'])
    else:
        return(r['PropertyTable']['Properties'][0]['CanonicalSMILES'])
def get_classyfire_data(dataset, classes= ['kingdom','superclass','class'],input_type='smiles'):
    kingdom=[]
    superclass=[]
    plainclass=[]
    for index, row in tqdm(dataset.iterrows(), total = len(dataset)):

        results_temp = GNPS_classyfire(row['reference_smiles'], classes=classes, input_type=input_type)
        kingdom.append(results_temp[0])
        superclass.append(results_temp[1])
        plainclass.append(results_temp[2])
    dataset['kingdom']=kingdom
    dataset['superclass']=superclass
    dataset['class']=plainclass
    return(dataset)
#


def GNPS_NPclassyfire(inputt, classes=None, input_type ='smiles'):
    # pass
    returning = []

    if classes ==None:
        print("no class type was fed")
        return()
    # if input_type =='inchikey':
    #     try:
    #         r = requests.get('https://gnps-classyfire.ucsd.edu/entities/%s.json'%inputt).json()
    #         if r != 'Key not found, try again later as we update our cache':
    #             for class_type in classes:
    #                 returning.append(r[class_type]['name'])
    #             return(returning)
    #         else:
    #             return("inchi_not_found")
    #     except:
    #         return("classyfire_error inchikey")
    if input_type == 'smiles':
        # print("i am in smiles if")
        # return(r)
        try:
            r = requests.get('https://npclassifier.ucsd.edu/classify?smiles=%s'%inputt).json()
        except:
            for class_type in classes:
                returning.append('classyfire_error_smiles')
            return(returning)
            # return("classyfire_error smiles")
        if r != {'message': 'unable to import structure'}:
            for class_type in classes:
                if r[class_type]!= None:
                    returning.append(r[class_type])
                else:
                    returning.append('class_not_found')

            return(returning)
                # return r[class_type]['name']
        else:
            for class_type in classes:
                returning.append('smiles_not_found')
            return(returning)
def GNPS_classyfire(inputt, classes=None, input_type ='smiles'):
    # pass
    returning = []
    
    if classes ==None:
        print("no class type was fed")
        return()
    # if input_type =='inchikey':
    #     try:
    #         r = requests.get('https://gnps-classyfire.ucsd.edu/entities/%s.json'%inputt).json()
    #         if r != 'Key not found, try again later as we update our cache':
    #             for class_type in classes:
    #                 returning.append(r[class_type]['name'])
    #             return(returning)
    #         else:
    #             return("inchi_not_found")
    #     except:
    #         return("classyfire_error inchikey")
    if input_type == 'smiles':
        # print("i am in smiles if")
        # return(r)
        try:
            r = requests.get('https://npclassifier.ucsd.edu/classify?smiles=%s'%inputt).json()
        except:
            for class_type in classes:
                returning.append('classyfire_error_smiles')
            return(returning)
            # return("classyfire_error smiles")
        if r != {'message': 'unable to import structure'}:
            for class_type in classes:
                if r[class_type]!= None:
                    returning.append(r[class_type]['name'])
                else:
                    returning.append('class_not_found')
                    
            return(returning)
                # return r[class_type]['name']
        else:
            for class_type in classes:
                returning.append('smiles_not_found')
            return(returning)
        



def GNPS_get(content, inputt='inchikey' , outputt='smiles',firstblock_only = False):
    from toolsets.std_list_prep import check_mol
    try:
        if firstblock_only == False:
            r = requests.get('https://gnps-structure.ucsd.edu/%s?%s=%s' % (outputt, inputt, content))
        elif inputt == 'inchikey' and firstblock_only == True:
            # print("i am in true block")/
            r = requests.get('https://gnps-structure.ucsd.edu/%s?%s=%s' % (outputt, inputt, content[0:14]))

        if r.text != '{"message":"structure cant be identified"}\n' and 'Server Error' not in r.text:
            return(r.text)
        else:
            # print("input %s not found" %inputt)
            return(np.NaN)
    except:

        print("something going wrong with API call, returning NaN")
        return(np.NaN)
    # print("i am in new method")
    # allowed values: inchikey, inchi, smiles
    # https://gnps-structure.ucsd.edu/smiles?inchi=<inchi string>
    # r = requests.get(f'https://gnps-structure.ucsd.edu/{outputt}?{inputt}={content}').json()


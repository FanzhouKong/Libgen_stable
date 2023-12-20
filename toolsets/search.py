import pandas as pd
import numpy as np
import numexpr
def string_search(data, column_name,item, reset_index = True,reverse = False):
    if reverse == False:
        _data= data[data[column_name].to_numpy() == item]
        # return data[data[column_name].to_numpy() == item]
    else:
        _data= data[data[column_name].to_numpy() != item]
    if reset_index == True:
        _data.reset_index(inplace= True, drop = True)
    return(_data)
        # return data[data[column_name].to_numpy() != item]
def quick_search_sorted(data_raw, column_name,value_start, value_end):
    # data.sort_values(by=column_name, inplace = True)
    search_array=data_raw[column_name].to_numpy(dtype="float")
    # data = data_raw.copy()
    index_start, index_end = search_array.searchsorted([value_start, value_end])
    # index_start = np.searchsorted(data[column_name], value_start,side = 'left')
    # index_end = np.searchsorted(data[column_name], value_end,side = 'right')
    return(data_raw.iloc[index_start:index_end])
def quick_search_values(data_raw, column_name,value_start, value_end, ifsorted = False):

    if ifsorted == False:
        data_sorted = data_raw.sort_values(by=column_name)
    data_return = quick_search_sorted(data_sorted, column_name, value_start, value_end)
    # index_start = np.searchsorted(data[column_name], value_start,side = 'left')
    # index_end = np.searchsorted(data[column_name], value_end,side = 'right')
    return(data_return)

def num_search(data, column_name,number, direction, step = None,inclusion = False):
    x = data[column_name].values
    if direction == ">":
        if inclusion == False:
            return(data[numexpr.evaluate('(x > number)')])
        else:
            return(data[numexpr.evaluate('(x >= number)')])
    elif direction == '<':
        if inclusion == False:
            return(data[numexpr.evaluate('(x < number)')])
        else:
            return(data[numexpr.evaluate('(x <= number)')])

    elif direction == '==':
        return(data[numexpr.evaluate('(x == number)')])
    elif direction =='between' and step != None:
        if inclusion == False:
            temp = data[numexpr.evaluate('(x > number-step)')]
            x = temp[column_name].values
            return (temp[numexpr.evaluate('(x < number+step)')])
        else:
            temp = data[numexpr.evaluate('(x >= number-step)')]
            x = temp[column_name].values
            return (temp[numexpr.evaluate('(x <= number+step)')])
    else:
        print('the wrong method is passed')




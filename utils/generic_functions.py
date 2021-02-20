import os
import pandas as pd


def file_exists(file_to_check):
    '''
    Input:
        file_to_check   # file name to check if exists
    Return:
        True if exists
    '''
    if os.path.isfile(file_to_check):
        return True
    return False

def directory_exists(dir_to_check):
    '''
    Input:
        dir_to_check    # directory to check
    Return:
        True if directory exists
    '''
    if os.path.isdir(dir_to_check):
        return True
    return False

def save_dataframe_to_excel_multiple_sheet(df, file, sheet_name):
    '''
    Input:
        df          # panda dataframe
        file        # file name to save dataframe
        sheet_name  # name of the sheet
    Return:
        True
    '''
    df.to_excel(file, sheet_name = sheet_name, index = False, header = True)
    return True

import os

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

# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 10:32:57 2024

@author: Collin
"""
# base_dir = User input file path to the base directory (where user stores other python scripts/data) 
# file_path = the specific file path to the dataset to be read in/create vars

def set_paths():
    # Prompt the user for the base directory path
    base_dir = input("Enter the base directory path: ")

    # Prompt the user for the specific file name
    file_name = input("Enter the data file's path: ")

    # Concatenate base_dir with the file name to create in_file path
    in_file = base_dir + file_name

    return base_dir, in_file
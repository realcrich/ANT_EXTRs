# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:26:32 2024

@author: Collin
"""
'''
from cytrack import *
from funcs.run_ARDetect import run_ARDetect

get_cytrack_main('C:/Users/Collin/anaconda3/Lib/site-packages/cytrack/cytrack_inputs.cfg')

run_ARDetect('D:/Research/DATA/wille_ar_detections/wille_ar_detections_inputs.txt')
'''
from cytrack import *
from funcs.run_ARDetect import run_ARDetect

# Run the cytrack main function
def run_cytrack_main():
    get_cytrack_main('C:/Users/Collin/anaconda3/Lib/site-packages/cytrack/cytrack_inputs.cfg')

# Run the ARDetect function
def run_ar_detect():
    run_ARDetect('D:/Research/DATA/wille_ar_detections/wille_ar_detections_inputs.txt')

if __name__ == "__main__":
    import multiprocessing

    p1 = multiprocessing.Process(target=run_cytrack_main)
    p2 = multiprocessing.Process(target=run_ar_detect)

    p1.start()
    p2.start()

    p1.join()
    p2.join()
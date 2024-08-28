from cytrack.cytrack_functions import get_dates, get_dates_vectors
from funcs.ARDetect_funcs import *

def run_ARDetect(args_file):
    
    """
    Run AR detection scheme from Wille et. al, 2019 and 
    plot polar projection maps of vIVT, cyclones, and diabatic heating 

    Parameters:
    args_file (dtype: str): name of file containing input parameters for run_ARDetect.py

    Returns:
    None (dtype: None): Runs wille_vivt_ar_detect.py in ARDetect_funcs.py using args_file
    """

    # Read arguments from the text file
    args = read_args_from_file(args_file)

    # Convert numeric values from strings to their appropriate types
    args['dt_h'] = int(args['dt_h'])
    args['prev_days'] = int(args['prev_days'])

    # Pass arguments to your functions
    idir= args['idir']
    dates, hours = get_dates_vectors(
        year_case_init=args['year_case_init'],
        month_case_init=args['month_case_init'],
        day_case_init=args['day_case_init'],
        hour_case_init=args['hour_case_init'],
        year_case_end=args['year_case_end'],
        month_case_end=args['month_case_end'],
        day_case_end=args['day_case_end'],
        hour_case_end=args['hour_case_end'],
        dt_h=args['dt_h'],
        prev_days=args['prev_days'],
        calendar=args['calendar']
        )
    plot_cyclones = args['plot_cyclones']
    plot_Q = args['plot_Q']
    
    wille_vivt_ar_detect(idir, dates, hours, plot_cyclones,plot_Q)
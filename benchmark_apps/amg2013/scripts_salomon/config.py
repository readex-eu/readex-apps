import collections
from collections import OrderedDict as ordDict

# Folder with all measurement data + DynDetect report
#
# Data type: List of tuples
#
# The first element in tuple is a path to the folder with measurement data,
# the second one is a path to an output TeX-file from DynDetect.
# If there's no DynDetect file to be included in the report, the second element
# is None.
#
# Examples:
#   root_folder = [('/home/user/callTreeEspreso', None)]
#   root_folder = [('/home/martin/Documents/READEX/measurements/TESTMPI/',
#                  '/home/martin/Documents/READEX/measurements/TESTMPI/report.tex')]
root_folder = [('../amg2013_meric_dir',None)]

# Main region rolea (used just as a title, no direct impact to visualized results)
# and name
#
# Data type: List of Dictionaries/Strings
#
# Key of the dictionary is a role, value is the region name itself. If no role is
# listed, then it's default role is "Main reg".
#
# Examples:
#   main_reg = [{'Solver': 'Solve_RegCG_singular_dom_itersolver.cpp:109'}]
#   main_reg = ['Solve_RegCG_singular_dom_itersolver.cpp:109']
main_reg = ['main_phase']

# Listing of visualized quantities/variables
#
# Data type: OrderedDictionary
#
# Key of dictionary is the category of visualized quantity
# (e.g. 'Blade summary), the value is a list of tuples,
# where the first element is a label of visualized q.
# (e.g. 'Energy consumption [J] (Samples)') and the second one
# is a unit of the quantity.
#
# Example:
#   y_label = ordDict((('Blade summary', [('Energy consumption [J] (Stats structure)', 'J')]),
#                     ('Job info - hdeem', [('Runtime of function [s]', 's')])))
y_label = ordDict((('COUNTERS - RAPL:', [('AVG Energy summary [J]', 'J')]),
                  ('Job info - RAPL', [('AVG Runtime of function [s]', 's')])))

# Time and energy variable
#
# Data type: Dictionary/None
#
# When used, than reports will contain information about
# run-time change for energetically-optimal settings
#
# Example:
#   1) time_energy_vars = None
#   2) time_energy_vars = {'time': ('Job info - hdeem', 'Runtime of function [s]'),
#                          'energy': ('Blade summary', 'Energy consumption [J] (Samples)')}
time_energy_vars = {'time': ('Job info - RAPL', 'AVG Runtime of function [s]'),
                    'energy': ('COUNTERS - RAPL:', 'AVG Energy summary [J]')}

# Parameters format from the name of CSV file
#
# Parameter ORDER MUST CORRESPOND to the order of
# parameters in the name of CSV file!
#
# Positions:
#   xLabel - input variable -> visualized on x-axis in plots
#   funcLabel - visualized as a single line in a plot
#   config - single data category, every config elements
#            combination forms one section in a generated report
#   key - additional parameter - there can be an arbitrary amount of them,
#         they're NOT visualized in plots or heat-maps
#
# Data type: List of lists of Strings
#
# Example:
#  1) file_name_args_tup = [['xLabel', 'Number of threads'],
#                        ['funcLabel', 'Frequency'],
#                        ['config', 'Preconditioner'],
#                        ['config', 'Uncore frequency']
#                        ]
#
#     For CSV file looking like this:
#        xLabel_funcLabel_config_config.csv (e.g. 12_2500000_PREC1_24.csv)
#
#  2) file_name_args_tup = [
#                        ['config', 'Node'],
#                        ['key', 'Threads'],
#                        ['xLabel', 'Frequency'],
#                        ['funcLabel', 'Uncore freq']
#                       ]
#
# The only compulsory position, which must be listed in CSV
# files' name is xLabel.
file_name_args_tup = [
                   ['config', 'Node'],
                   ['key', 'threads'],
                   ['xLabel', 'CF'],
                   ['funcLabel', 'UCF']
                  ]

# Default values for keys' values
#
# Data type: List
#
# This value really MUST be used in some CSV file names
# (i.e. you can't have a default key value 1, when you've
# only measured with 5 and 10 threads etc.)
#
# Example:
#   1) def_keys_vals = []
#   2) def_keys_vals = ['24', '12']
def_keys_vals = ['12']

# Units for keys
#
# Data type: List
#
# Units of keys' values (see above)
#
# Example:
#   1) keys_units = []
#   2) keys_units = ['th', 'Hz']
keys_units = ['threads']

# Default value, unit and multiplier for xLabel (x-axis) values
#
# Data types: Integer, String, Float/None
#
# The default value really MUST be used in some CSV file names
# (i.e. you can't have a default value 24, when you've
# only measured with 12 and 16 threads etc.)
#
# The multiplier can be used to convert default value to
# other units.
#
# Example:
#   def_x_val = 25
#   x_val_unit = 'GHz (core)'
#   x_val_multiplier = 0.1
def_x_val = 25
x_val_unit = 'GHz CF'
x_val_multiplier = 0.1

# Default value, unit and multiplier for funcLabel
# (single lines in plots) values
#
# Data types: Integer, String, Float/None
#
# The default value really MUST be used in some CSV file names
# (i.e. you can't have a default value 2500000, when
# you've only measured with 1200000 and 1600000 Hz etc.)
#
# The multiplier can be used to convert default value to
# other units.
#
# Example:
#   def_label_val = 30
#   func_label_unit = 'GHz (uncore)'
#   label_val_multiplier = 0.1
def_label_val = 30
func_label_unit = 'GHz UCF'
label_val_multiplier = 0.1

# Nested functions
#
# Data type: Ordered Dictionary of tuples in format (String, list of strings) /
#            List of strings
#
# Ordered dictionary of tuples, where the first element of tuple
# is the name of position (e.g. Preconditioner) and the second
# element is the list of region names on that position
# (e.g. ['Projector_l_inv_compG_itersolver.cpp:386',
#        'Projector_l_inv_compG_itersolver.cpp:400']).
#
# If you don't want to list positions specifically and you're ok with positions
# named the same as regions itself, it's enough to provide just list with regions' names.
#
# Examples:
#  1) allNestedRegsDic = collections.OrderedDict((('F operator', ['apply_A_l_comp_dom_B_itersolver.cpp:446']),
#                                                 ('Preconditioner', ['apply_prec_comp_dom_B_itersolver.cpp:394']),
#                                                 ('Projector', ['Projector_l_inv_compG_itersolver.cpp:386',
#                                                                'Projector_l_inv_compG_itersolver.cpp:400'])))
#
#  2) all_nested_regs = ['B']
all_nested_regs = ['set_up','hypre_pcgsetup','hypre_pcgsolve']

smooth_runs_average = None

# Name of the iteration region
#
# This can be the region nested in the main one, but superior
# to the nested ones. It's calls are visualized as single
#  iterations.
#
# If you don't want to specify it, assign None value (dynamic
# behavior will be evaluated according to single main region
# calls).
#
# Data type: String / None
#
# Example:
#   iter_call_region = 'iter'
iter_call_region = None

# Turn on/off generating report part
# with single program starts
#
# Data type: Boolean
#
# Example:
#   detailed_info = True
detailed_info = True

# Path to the .options file with individual settings
# for specific regions
#
# Data type: String / None
#
# Example:
#   options_file_path = None
options_file_path = None

optim_settings_file_path = None

# Dictionary of frequency and uncore frequency
# variables according to the 'file_name_args_tup'
# (see above)
#
# Turn on/off generation of .opts files
# with the optimal settings for every
# nested region
#
# Format: {'value_from_file_name_args_tup':'FREQUENCY',
#          'other_value_from_file_name_args_tup':'UNCORE_FREQUENCY'}
#
# Data type: Dictionary/None
#
# Examples:
#   1) generate_optim_settings_file = None
#   2) generate_optim_settings_file = {'Frequency': 'FREQUENCY',
#                                      'Uncore freq': 'UNCORE_FREQUENCY'}
generate_optim_settings_file = {'threads': 'NUM_THREADS',
                             'CF': 'FREQUENCY',
                             'UCF': 'UNCORE_FREQUENCY'}

# RAPL baseline energy consumption
# (eng per time, e.g. 100J per second)
#
# Data type: Integer / String / None
#
# WARNING: If baseline is an arithmetic
# formula, then used configuration,
# keys, func label and X label
# must be numbers!
#
# 'config.time_energy_vars' must be set
# when this variable is not None!
#
# Arithmetic formula:
#   conf[i] - i-th configuration element
#   keys[i] - i-th key
#   func_lab - function label
#   x_lab - X label
#
# Example:
#   baseline = 100
#   baseline = 'func_lab * 2 + keys[0]'
baseline = 70

test_csv_init = False

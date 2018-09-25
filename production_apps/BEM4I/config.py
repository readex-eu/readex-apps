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
#   rootFolder = [('/home/user/callTreeEspreso', None)]
#   rootFolder = [('/home/martin/Documents/READEX/measurements/TESTMPI/',
#                  '/home/martin/Documents/READEX/measurements/TESTMPI/report.tex')]

#root_folder = [('/home/zapletal/bem4i_readex/bem4i_meric',None)]
root_folder = [('./bem4i_meric_dir',None)]

# Main region rolea (used just as a title, no direct impact to visualized results)
# and name
#
# Data type: List of Dictionaries/Strings
#
# Key of the dictionary is a role, value is the region name itself. If no role is
# listed, then it's default role is "Main reg".
#
# Examples:
#   mainReg = [{'Solver': 'Solve_RegCG_singular_dom_itersolver.cpp:109'}]
#   mainReg = ['Solve_RegCG_singular_dom_itersolver.cpp:109']

main_reg = ['main_phase']

# Listing of visualized quantities
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
# yLabel = ordDict((('Blade summary', [('Energy consumption [J] (Stats structure)', 'J')]),
#                  ('Job info - hdeem', [('Runtime of function [s]', 's')])))

y_label = ordDict((('Blade summary', [('Energy consumption [J] (Samples)', 'J')]),
                  ('Job info - hdeem', [('Runtime of function [s]', 's')])))

time_energy_vars = {'time': ('Job info - hdeem', 'Runtime of function [s]'),
                    'energy': ('Blade summary', 'Energy consumption [J] (Samples)')} 

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
#  1) fileNameArgsTup = [['xLabel', 'Number of threads'],
#                        ['funcLabel', 'Frequency'],
#                        ['config', 'Preconditioner'],
#                        ['config', 'Uncore frequency']
#                        ]
#
#     For CSV file looking like this:
#        xLabel_funcLabel_config_config.csv (e.g. 12_2500000_PREC1_24.csv)
#
#  2) fileNameArgsTup = [
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
# Data type: List of Strings
#
# This value really MUST be used in some CSV file names
# (i.e. you can't have a default key value 1, when you've
# only measured with 5 and 10 threads etc.)
def_keys_vals = ['24']

# Units for keys
#
# Data type: List of Strings
keys_units = ['threads']

# Default value for xLabel (x-axis) values
#
# MANDATORY PARAMETER!
#
# Data type: Integer
#
# This value really MUST be used in some CSV file names
# (i.e. you can't have a default value 24, when you've
# only measured with 12 and 16 threads etc.)
def_x_val = 25

# Unit for x-value
#
# MANDATORY PARAMETER!
#
# Data type: String
x_val_unit = 'GHz CF'

# Multiplier for x-value
#
# Data type: Float / None
x_val_multiplier = 0.1

# Default value for funcLabel (single lines in plots)
# values
#
# Data type: Integer / None
#
# This value really MUST be used in some CSV file names
# (i.e. you can't have a default value 2500000, when
# you've only measured with 1200000 and 1600000 Hz etc.)
def_label_val = 30

# Units for funcLabel values
#
# Data type: String / None
func_label_unit = 'GHz UCF'

# Label-val multiplier
#
# Data type: Float / None
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
#  2) allNestedRegs = ['B']

all_nested_regs = ['assemble_v','assemble_k','gmres_solve','print_vtu']

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
#   iterCallRegion = 'iter'
#   iterCallRegion = None

iter_call_region = None
#iterCallRegion = 'gmres_iter'

# Turn on/off generating report part
# with single program starts
#
# Data type: Boolean
#
# Example:
#   detailedInfo = True

detailed_info = True
#detailedInfo = False

options_file_path=None

# Variable which turns on/off generating
# of .opts file with the optimal settings
# for every region.
#
# Dictionary containing "mapping"
# from fileNameArgsTup values (NOT keys)
# to:
#   FREQUENCY
#   UNCORE_FREQUENCY
#   NUM_THREADS
#
# Or None if you don't want to generate
# the .opts file.
#
# Data type: Dictionary / None
#
# Example:
#   generateOptimSettingsFile = {'Threads': 'NUM_THREADS',
#                                'CF': 'FREQUENCY',
#                                'UCF': 'UNCORE_FREQUENCY'}
#
#   generateOptimSettingsFile = None

generate_optim_settings_file = {'threads': 'NUM_THREADS',
                             'CF': 'FREQUENCY',
                             'UCF': 'UNCORE_FREQUENCY'}

# Path to the .opts file with individual settings
# for specific regions
#
# If you want to use this variable, you HAVE TO
# set the dictionary in 'generateOptimSettingsFile'!
#
# Data type: String/None
optim_settings_file_path = None

baseline=None

#
# script used to run the mass fit in the different neutron classes, for a given binning and process
#

import os
import getopt, sys

# Remove 1st argument from the
# list of command line arguments
argumentList = sys.argv[1:]

#input options
options = "hb:i:"

#input long options
long_options = ["help","bins=","id="]

try:
  #default values
  phiBins = '1'
  identifier = 'jPsi'

  # Parsing argument
  arguments, values = getopt.getopt(argumentList, options, long_options)
  # checking each argument
  for currentArgument, currentValue in arguments:
    if currentArgument in ("-h","--help"):
      print('----------------------------')
      print('Usage: python3 run-mass-fits.py --bins <phi_bins> --id <process_name>')
      print('Available options:')
      print('--bins    number of phi bins')
      print('--id      id of the process')
      print('----------------------------')
      sys.exit()
    
    elif currentArgument in ("-b","--bins"):
      phiBins = str(currentValue)

    elif currentArgument in ("-i","--id"):
      identifier = str(currentValue)

  nClasses = ["0n0n","Xn0n","XnXn"]
  for nc in nClasses:
    comm = 'root -l -q \'fitJPsiInPhiBins.C+("' + str(nc) + '",' + str(phiBins) + ',"' +  str(identifier) + '")\''
    print(comm)
    os.system(comm)

except getopt.error as err:
  print("error: " + str(err))
  print("use [-h, --help] to see available options")
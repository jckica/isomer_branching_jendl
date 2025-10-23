"""
******************************************************************************
  Script used to generate a depletion chain
******************************************************************************
"""

#******************************************************************************
# PACKAGES
#******************************************************************************

import openmc
import openmc.deplete
import os
import glob
from datetime import datetime

#******************************************************************************
# FUNCTIONS
#******************************************************************************

#- Function for printing the date
def print_time(message):
    print()
    print(message+str(datetime.now()))

#******************************************************************************
# GENERATION OF DEPLETION CHAIN
#******************************************************************************

print_time('===> Running started on ')

#- Load ENDF sub-libraies: decay, fission product yield & neutrons
decay_files = glob.glob('./JENDL/jendl5-dec_upd5/jendl5-dec_upd5/*.dat')
print(decay_files)
fpy_files = glob.glob('./JENDL/jendl5-fpy_upd8/jendl5-fpy_upd8/*.dat')
print(fpy_files)
neutron_files = glob.glob('./JENDL/jendl5-n/jendl5-n/*.dat')
print(neutron_files)

#- Generate depletion chain
chain = openmc.deplete.Chain.from_endf(decay_files,fpy_files,neutron_files,reactions=list(openmc.deplete.chain.REACTIONS.keys()),progress=True)
chain.export_to_xml('JENDL_chain.xml')

print_time('===> Processing completed on ')
print()
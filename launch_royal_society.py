#!/usr/bin/python
from os import system,popen
import random
import sys
import time

############################################
## SYSTEMATICAL TEST OF ARIDITY/MORTALITY ##
############################################

# Two parameters are varied:
# 1) the mortality induced by water limitation (called "water_mr" here)
# 2) the soil aridity (somewhat confusingly called "water_needed" here)

## in the final version of the paper, we tested the following ranges:
# 1) i in [0,1], with a step of 0.1, for the mortality (water_mr)
#    With the argument parameter recoding, this is achieved with: i_range = range(11)
# 2) j in [0,1], with a step of 0.1, for the soil aridity (water_needed)
#    With the argument parameter recoding, this is achieved with: j_range = range(1,21,2)

i_range = range(11)
j_range = range(1,21,2)

#################
# Change here to launch the simulations with different mortality/growth mechanisms:
#################

config = "watergrid40h_2m.ini" 
outputprefix = "mediumsmall_grid40h_2m_"

#config = "watergrid40h_21.ini" 
#outputprefix = "mediumsmall_grid40h_21_"

#config = "watergrid40h_mm.ini" 
#outputprefix = "mediumsmall_grid40h_mm_"

#config = "watergrid40h_m1.ini" 
#outputprefix = "mediumsmall_grid40h_m1_"



argc = len(sys.argv)
if argc > 1:
  extra_args = ' '.join(sys.argv[1:argc])
else:
  extra_args = ''

exe = "mod"

random.shuffle(i_range)
random.shuffle(j_range)
for i in i_range:
  for j in j_range:
      time.sleep(1)
      # Wait for available computational resources:
      while (float(popen("mpstat 1 1 | grep 'Average' | awk '{ printf $NF }'").read()) < 5):
          print "==>", (str(float(popen("mpstat 1 1 | grep 'Average' | awk '{ printf $NF }'").read()))), "<=="
          time.sleep(35)
      print "will start a new process", "==>", (str(float(popen("mpstat 1 1 | grep 'Average' | awk '{ printf $NF }'").read()))), "<=="
      # argument recoding:
      water_needed = str(j/10.0*3.0) # scaled in [0,6] in C++
      water_mr = str(i/10.0) # scaled in [0,1] in C++
      outp = outputprefix + water_mr + "_" + water_needed # intermediate steps
      system("screen -dmS 'ibm"+str(j)+"' bash -c 'LD_PRELOAD=./libpngwriter.so ./"+exe+" -c "+config+" --WATER_AVAILABLE "+water_needed+" --S0.water_mr "+water_mr+" --S1.water_mr "+water_mr+" --output "+outp+" "+extra_args+"'") 



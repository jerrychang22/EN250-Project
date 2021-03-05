#!/usr/bin/python
from os import system,popen
import random
import sys
import time

###############################
## SYSTEMATICAL TEST OF ARIDITY
###############################

# the soil aridity ("water_needed") is varied, while the mortality is kept constant

i_range = [2]
j_range = range(1,21,1)

#################
# Change here to launch the simulations with implicit and explicit water modeling:
#################

config = "water_explicit.ini"
#config = "water_implicit.ini"
#config = "water_hires_explicit.ini"
#config = "water_hires_implicit.ini"

outputprefix = "water_explicit_bis_"
#outputprefix = "water_implicit_bis_"
#outputprefix = "water_1_hires"
#outputprefix = "water_2_hires"



argc = len(sys.argv)
if argc > 1:
  extra_args = ' '.join(sys.argv[1:argc])
else:
  extra_args = ''

exe = "mod2"

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
      rnd = str(random.randint(1,10**10))
      print("screen -dmS 'ibm"+str(j)+"_"+str(rnd)+"' bash -c 'LD_PRELOAD=./libpngwriter.so ./"+exe+" -c "+config+" --WATER_AVAILABLE "+water_needed+" --S0.water_mr "+water_mr+" --S1.water_mr "+water_mr+" --output "+outp+" "+extra_args+"'") 
      system("screen -dmS 'ibm"+str(j)+"_"+str(rnd)+"' bash -c 'LD_PRELOAD=./libpngwriter.so ./"+exe+" -c "+config+" --WATER_AVAILABLE "+water_needed+" --S0.water_mr "+water_mr+" --S1.water_mr "+water_mr+" --output "+outp+" "+extra_args+"'") 


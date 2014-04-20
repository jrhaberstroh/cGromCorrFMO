import chimera
from subprocess import call
import sys

#print sys.argv
#
call(["python","f2_2GenerateModeColormap.py",sys.argv[1],sys.argv[2],'modes.com'])

chim = chimera.openModels.open("../../data/newFMO.gro")

#with open('../../data/chimera/FMO_1_1_chimera.com','r') as f:
with open('modes.com','r') as f:
    for l in f:
        print l
        chimera.runCommand(l.strip())


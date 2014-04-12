import MDAnalysis

home = '/home/jhaberstroh'
top = home + "/Code/photosynth/FMO_Forcefield/scripts/mddir/4BCL.top"
gro = home + "/Code/photosynth/FMO_Forcefield/scripts/mddir/4BCL.gro"

uni = MDAnalysis.Universe(top, gro)

import MDAnalysis

mddir = "/home/jhaberstroh/mass-storage/data/2014-08-4BCL/md/md_0"
md_data = MDAnalysis.Universe(mddir+"/md_100ps.tpr", mddir+"/md_100ps.trr")
BCL_atoms = md_data.selectAtoms('byres (resname BCL)')

for i,ts in enumerate(md_data.trajectory):
    if i > 10:
        break
    bcl_pos = np.zeros((len(bcl_atoms), 3))
    for j,atom in enumerate(bcl_atoms):
        bcl_pos[j,:] = atom.pos[:]
        print atom.pos

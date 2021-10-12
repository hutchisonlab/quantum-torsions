#!/usr/bin/env python

# Thanks to Peter Schmidtke for the guts of this script
# https://pschmidtke.github.io/blog/rdkit/crystallography/small%20molecule%20xray/xray/database/2021/01/25/cod-and-torsion-angles.html

import pandas as pd
import numpy as np
import glob

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

# ignore warnings
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')

# from ETKDG paper: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00654
torsions=pd.read_table("ring_smarts_patterns.txt",delim_whitespace=True,header=None,usecols=[1])
# from Ring ETKDG paper: https://pubs.acs.org/doi/10.1021/acs.jcim.0c00025
# torsions=pd.read_table("ring_smarts_patterns.txt",header=None,usecols=[1])

# filename template for output, e.g. t1.txt, t2.txt, etc.
out_template = 't{}-rings.txt'

# print 
# patterns=torsions[1]
# for debugging
patterns=torsions[1][:-1]

# bin size (in degrees)
bin_size = 1
bins = round(360 / bin_size)

# This outer loop is for each one of the torsion patterns
# We compile the pattern, then loop through the SDF / XYZ files
index = 0
for torsionSmarts in patterns:
    index += 1
    print(index, torsionSmarts)

    angles = np.zeros(bins) # create a histogram of angles with X degree bins
    total_matches = 0

    torsionQuery = Chem.MolFromSmarts(torsionSmarts)

    # these SMARTS have atom maps, so convert them
    # http://www.rdkit.org/docs/GettingStartedInPython.html#atom-map-indices-in-smarts
    index_map = {}
    for atom in torsionQuery.GetAtoms() :
        map_num = atom.GetAtomMapNum()
        if map_num:
            index_map[map_num-1] = atom.GetIdx()

    map_list = [index_map[x] for x in sorted(index_map)]
    
    try:
        # loop through the files and then each of the molecules in the file (if more than one)
        for sd_file in glob.iglob("/zfs1/ghutchison/dlf57/crest-zinc/*/*-conf.sdf"):
            suppl = Chem.SDMolSupplier(sd_file, removeHs=False)

            for mol in suppl:

                if mol is None: continue

                # get the 3D geometry
                conf = mol.GetConformer(0)

                matches = mol.GetSubstructMatches(torsionQuery)
                for match in matches:
                    total_matches += 1 # to normalize

                    # get the atom maps from the SMARTS match
                    mapped = [match[x] for x in map_list]
                    angle = rdMolTransforms.GetDihedralDeg(conf, mapped[0],mapped[1],mapped[2],mapped[3])
                    if (angle < 0.0):
                        angle += 360.0

                    # okay, we want to hash - e.g., 5° bins    
                    angle = round(angle / bin_size) % bins

                    angles[angle] += 1

        if total_matches > 0:
            # angles = angles / total_matches # normalize
            np.savetxt(out_template.format(index), angles, fmt='%.3e', delimiter=',')
            # print(angles)
    except ValueError:
        print('Issue with {} - {}'.format(index, torsionSmarts))
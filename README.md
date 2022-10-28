# Quantum Torsions

Supporting Data for "Systematic Comparison of Experimental Crystallographic Geometries and Gas-Phase Computed Conformers for Torsion Preferences"

by Dakota Folmsbee (@dlf57) David Koes (@dkoes) and Geoffrey Hutchison (@ghutchis)

All CREST / GFN2 low-energy conformations (over 2.5GB of compressed `.tar.bz2` files) are available through Figshare:
https://doi.org/10.6084/m9.figshare.21395061.v1

Uncompressed, this data can occupy over 20GB of disk space. GitHub does not permit large files.

## tldr: Use qtdg.ipynb and qtdg.txt to sample torsions

We have included [`qtdg.txt`](https://github.com/hutchisonlab/quantum-torsions/blob/main/qtdg.txt) as acyclic torsions SMARTS patterns and [`qtdg.ipynb`](https://github.com/hutchisonlab/quantum-torsions/blob/main/qtdg.ipynb)
as an example RDKit notebook to do torsion sampling / driving using the 
derived quantum torsion distributions from this work.

In short, generate a uniform random number, and you'll get back a dihedral angle in degrees from 0..360.

These are generated from `gauss-fit.txt` using [`inverse-data.ipynb`](https://github.com/hutchisonlab/quantum-torsions/blob/main/inverse-data.ipynb) which calculates the cumulative probability distribution from Gaussian fits, normalizes and inverts to generate a set of linear interpolations for dihedral
angles (in degrees from 0..360) from uniform random numbers.

## Inventory of Files and Directories

* `combined-histograms.ipynb`
* `combined-ring-histograms.ipynb`

Jupyter notebooks to generate fits, including `fit-figures`, comparisons with ETKDG etc. (e.g., `cos-fit.txt`, `cos-ring-fit.txt`, `gauss-fit.txt`, `gauss-ring-fit.txt`, etc.)

We include a cosine-squared fit which looked better on some patterns but was ultimately worse.

* `figures.ipynb`
* `ring-figures.ipynb`

Notebooks to generate `figures` and `ring-figures` including histograms of COD, CREST conformers for COD, and combined sets. Uses images from `smarts-figures` which are generated by https://smarts.plus/

* `cod-torsions`
Text files contain percent of torsions at that degree.
* `pqr-torsions`, `pubchemqc-torsions`, and `zinc-torsions`
Text files contain the raw number of torsions at that degree.
* `*-crest.tar.bz2`
Raw CREST output files in XYZ, usually with associated SDF/mol file

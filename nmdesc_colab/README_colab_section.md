## Run on Colab

The full analysis runs in a browser with no local setup:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CobanAkdemirlab/NMDescapediseasegene_paper/blob/main/NMDesc_colab.ipynb)

The notebook installs the R dependencies from a dated CRAN snapshot, clones this
repository, downloads the input files from a public archive, and runs the v4 pipeline.
Set `DATA_URL` in the first code cell to the archive location.

Inputs are not tracked in git: the v4 scripts read 352 distinct files, all derived from
public sources (ClinVar, gnomAD, Ensembl, OMIM). `colab_input_inventory.csv` lists them,
and `make_data_bundle.sh` assembles the archive from a local copy.

# GGpy

GGI automatization and feature extraction

Software requierements:

* pip
* python3


## Installation

It is advisable to install this package inside a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [python environment](https://docs.python.org/3/library/venv.html) if there are issues with system permissions.

Using `pip`:

<!-- # pip install numpy # needed for python<3.7 -->

```Bash
pip install ggi
```

Using `git` and `pip` (Optional):
<!-- python3 -m pip install numpy # needed for python<3.7 -->
```Bash
git clone https://github.com/Ulises-Rosas/GGpy.git
cd GGpy
python3 -m pip install .
```

## Usage

Main Command:

```Bash
ggpy -h
```

```
usage: ggpy [-h] {ggi,features,post} ...

                                 GGI and more
                                      

positional arguments:
  {ggi,features,post}
    ggi                Gene-Genealogy Interrogation (GGI)
    features           Features from both alignment and tree information

optional arguments:
  -h, --help           show this help message and exit
```
### GGI

```Bash
ggpy ggi demo/*fasta -t demo/ggi_tax_file.csv -H demo/myhypothesis.trees  
cat out_ggi.txt
```
```
alignment	tree_id	group	rank	au_test
E0055.fasta	1	(Outgroup,(Eso_salmo,(Argentiniformes,(Osme_Stomia,(Galaxiiformes,Neoteleostei)))));	1	0.880
E0055.fasta	2	(Outgroup,((Eso_salmo,Argentiniformes),(Osme_Stomia,(Galaxiiformes,Neoteleostei))));	2	0.120
E1532.fasta	2	(Outgroup,((Eso_salmo,Argentiniformes),(Osme_Stomia,(Galaxiiformes,Neoteleostei))));	1	0.921
E1532.fasta	1	(Outgroup,(Eso_salmo,(Argentiniformes,(Osme_Stomia,(Galaxiiformes,Neoteleostei)))));	2	0.079
```

Utilities

* `root_groups.py` : Root groups at ggpy results

### Features

Feature extraction from alignments and sequences

```Bash
ggpy features -A [alignment file extension] -T [tree file extension]
```
### TODO

- Add a machine learning model to make a non-linear regression between labels from ggi and features


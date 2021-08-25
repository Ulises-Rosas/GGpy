# GGpy

GGI automatization

Software requierements:

* pip
* python3


## Installation

Using `pip`:

```Bash
pip install numpy # needed for python<3.7
pip install ggi
```


Using `git` and `pip` (Optional):
```Bash
git clone https://github.com/Ulises-Rosas/GGpy.git
cd GGpy
python3 -m pip install numpy # needed for python<3.7
python3 -m pip install .
```

## Usage

Main Command:

```Bash
ggpy -h
```

```
usage: ggpy [-h] {ggi,features,post} ...

                                 GGI and Post-GGI
                                      

positional arguments:
  {ggi,features,post}
    ggi                Gene-Genealogy Interrogation (GGI)
    features           Features from both alignment and tree information
    post               Classification of GGI hypothesis based on features

optional arguments:
  -h, --help           show this help message and exit
```
### GGI

```Bash
ggpy ggi demo/*fasta -t demo/ggi_tax_file.csv -H demo/myhypothesis.trees  
cat out_ggi.txt
```
```
alignment       tree_id group   rank    au_test
E0055.fasta     1       (Outgroup,(Eso_salmo,(Argentiniformes,(Osme_Stomia,(Galaxiiformes,Neoteleostei)))));    1       0.882
E0055.fasta     2       (Outgroup,((Eso_salmo,Argentiniformes),(Osme_Stomia,(Galaxiiformes,Neoteleostei))));    2       0.118
E1532.fasta     2       (Outgroup,((Eso_salmo,Argentiniformes),(Osme_Stomia,(Galaxiiformes,Neoteleostei))));    1       0.922
E1532.fasta     1       (Outgroup,(Eso_salmo,(Argentiniformes,(Osme_Stomia,(Galaxiiformes,Neoteleostei)))));    2       0.078
```

Utilities

* `root_groups.py` : Root groups at ggpy results

### Features

```Bash
ggpy features -A [alignment file extension] -T [tree file extension]
```

### post-GGI

```Bash
ggpy post [ggi result] -f [features result] -c [comparison file]
```

# GGpy: 

GGI automatization

Software requierements:

* pip

## Installation

Using `pip`:

```Bash
pip install ggi
```

Development version by using `git` and `pip` (optional):
```Bash
git clone https://github.com/Ulises-Rosas/GGpy.git
cd GGpy
python3 -m pip install .
```

## Usage

```Bash
ggpy demo/*fasta -t demo/ggi_tax_file.csv -H demo/myhypothesis.trees  
cat out_ggi.txt
```
```
alignment       tree_id group   rank    au_test
E0055.fasta     1       (Outgroup,(Eso_salmo,(Argentiniformes,(Osme_Stomia,(Galaxiiformes,Neoteleostei)))));    1       0.882
E0055.fasta     2       (Outgroup,((Eso_salmo,Argentiniformes),(Osme_Stomia,(Galaxiiformes,Neoteleostei))));    2       0.118
E1532.fasta     2       (Outgroup,((Eso_salmo,Argentiniformes),(Osme_Stomia,(Galaxiiformes,Neoteleostei))));    1       0.922
E1532.fasta     1       (Outgroup,(Eso_salmo,(Argentiniformes,(Osme_Stomia,(Galaxiiformes,Neoteleostei)))));    2       0.078
```

```Bash
ggpy -h
```

```
usage: ggpy [-h] [-H] [-e] [-w] [-r] [-n] [-s] [-E] [-c] [-i] -t
            [alns [alns ...]]

                        Gene-Genealogy Interrogation

    * Standard usage:

        $ ggpy [exon files] -t [taxonomy file]

            note 1: The taxnomy file is CSV-formated and must 
                    contain the following:

                    names               group              
                    [sequence header 1],[group of header 1]
                    [sequence header 2],[group of header 2]
                    [sequence header 3],[group of header 3]
                     ...                 ...      

            note 2: For n groups in the taxonomy file all possible
                    topologies are generated if a topology file is 
                    not given (see below). Here is a table 
                    for reference

                    n | Unrooted trees | Rooted tree 
                    --+----------------+------------ 
                    3 |              - |           3 
                    4 |              3 |          15 
                    5 |             15 |         105 
                    6 |            105 |         945 
                    7 |            945 |      10 395 
                    8 |         10 395 |     135 135 
                    9 |        135 135 |           - 

    * Specify pre-defined hypothesis:

        $ ggpy [exon files] -t [taxonomy file] -H [file with topologies]
            
            note: `-h` accepts a group-based topologies

    * Specify pre-defined extended hypothesis:

        $ ggpy [exon files] -t [taxonomy file] -H [file with topologies] -e

            note: `-e` specifies that topologies at `-h` files are 
                  extended (i.e., with actual species and not groups)

required arguments:
  alns                  Alignments
  -t , --tax_file       Taxonomy file. Format in csv: "[sequence
                        name],[group]"

RAxML constrained tree parameters:
  -E , --evomol         [Optional] RAxML evol. model for constrained tree
                        inference [Default: GTRGAMMA]
  -c, --codon_aware     [Optional] If selected, codon partition file is added
  -i , --iterations     [Optional] Number of iterations for MLEs [Default: 1]

optional arguments:
  -h, --help            show this help message and exit
  -H , --hypothesis     Pre-defined hypothesis file, either extended or not
  -e, --extended        [Optional] If selected, file with topologies contains
                        extended trees (i.e., species instead of groups)
  -w, --write_extended  [Optional] If selected, extended topologies are
                        written and exit
  -r, --rooted          [Optional] If selected, all posible rooted topologies
                        are generated when pre-defined hypothesis file is not
                        given
  -n , --threads        [Optional] number of cpus [Default = 1]
  -s , --suffix         [Optional] Suffix for each written file [Default:
                        ggi.txt]
```

Utilities

* `root_groups.py` : Root groups at ggpy results

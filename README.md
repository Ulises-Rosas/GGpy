# GGpy

GGI automatization and feature extraction

Software requierements:

* pip
* python3


## Installation

It is advisable to install this package inside a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [python environment](https://docs.python.org/3/library/venv.html) if there are issues with system permissions.

### Option 1: Using `pip`:


```Bash
pip install numpy # needed for those with python<3.7
pip install ggi
```

### Option 2: Using `git` and `pip`:

```Bash
git clone https://github.com/Ulises-Rosas/GGpy.git
cd GGpy
python3 -m pip install numpy # needed for those with python<3.7 
python3 -m pip install .
```

## Usage

Main Command:

```Bash
ggpy -h
```

```
usage: ggpy [-h] {ggi,features,fi} ...

                                    GGI and more
                                      

positional arguments:
  {ggi,features,fi}
    ggi              Gene-Genealogy Interrogation
    features         Features from both alignment and tree information
    fi               Feature importance using Random Forest

optional arguments:
  -h, --help         show this help message and exit
```
### GGI

```Bash
ggpy ggi demo/*fasta -t demo/ggi_tax_file.csv -H demo/myhypothesis.trees -c
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

Feature extraction from alignments and sequences. List of features can be seen [here](https://github.com/Ulises-Rosas/fishlifetraits/blob/main/fishlifetraits/var_desc.md).


```Bash
ggpy features -A '.fasta' -T '.tree' --path demo/alns_trees
cat features_stats.tsv
```
```
aln_base  nheaders  pis vars  seq_len seq_len_nogap nogap_prop  gc_mean gc_var  gap_mean  gap_var pi_mean pi_std  total_tree_letreeness inter_len_mean  inter_len_var ter_len_mean  ter_len_var supp_mean coeffVar_len  rcv treeness_o_rcv  saturation  LB_std  invariants  singletons  patterns  entropy gc_mean_pos1  gc_var_pos1 gap_mean_pos1 gap_var_pos1  gc_mean_pos2  gc_var_pos2 gap_mean_pos2gap_var_pos2 gc_mean_pos3  gc_var_pos3 gap_mean_pos3 gap_var_pos3
E0011.fasta 63  85  96  183 111 0.981785  62.04049  13.83979  1.821494  22.029381 78.938945 7.879053  4.208196  0.578473  0.040572  0.001893  0.028157  0.001832  69.736667 1.269078  0.072226  0.124856  1.539796  16.08387  87  11  111 0.456799  69.191997 7.767883  1.821494  22.029381 39.702454 3.482264  1.821494  22.029381 77.227019 128.769901  1.821494  22.029381
E0010.fasta 27  67  99  177 80  0.919858  60.068892 33.133186 8.014229  199.349264  72.86365  17.616112 2.35872 0.517771  0.050887  0.010657  0.042127  0.00315 52.920833 1.751212  0.142243  0.274722  0.904961  36.806211 78  32  124 0.36323 62.862803 39.906434 8.474576  194.020286  34.739446 4.653384  7.784055  202.253831  82.62254  134.257174  7.784055  202.253831
E0001.fasta 76  93  115 189 120 0.991437  58.931207 9.757202  0.856307  13.758235 85.068412 5.965092  5.426492  0.444095  0.033012  0.002317  0.040765  0.005957  63.165278 1.742215  0.066273  0.149232  1.92  28.254218 74  22  131 0.345902  62.618231 7.868368  0.856307  13.735839 30.497204 4.321067  0.877193  13.867119 83.665039 56.482223 0.835422  13.738049
E0003.fasta 62  100 113 186 63  0.967135  56.586615 8.632091  3.286507  73.43884  74.963182 10.877145 4.4636  0.464956  0.035176  0.002313  0.03852 0.004174  68.249153 1.543607  0.092483  0.198907  0.812799  28.17858  73  13  135 0.503376  65.372518 5.371161  3.225806  73.096672 35.326784 6.965671  3.329865  73.256254 69.048492 72.207852 3.30385 74.156652
E0009.fasta 47  88  109 183 93  0.978723  53.971619 17.914953 2.12766 34.784089 81.438458 9.201035  3.226205  0.453635  0.033262  0.002115  0.037504  0.002676  60.656818 1.376859  0.081267  0.179146  1.118977  33.739651 74  21  127 0.388955  52.532865 9.785306  2.12766 34.926901 28.327448 4.880208  2.023021  35.078551 81.119742 129.053417  2.232299  34.63603
E0013.fasta 86  102 110 186 45  0.969242  62.607592 17.994162 3.075769  76.857924 78.725564 10.834178 4.610766  0.470344  0.026128  0.001177  0.028397  0.002364  64.292771 1.542837  0.09687 0.205956  0.607839  22.414153 76  8 144 0.416121  63.810444 11.495023 3.075769  76.857924 49.644551 5.542949  3.075769  76.857924 74.367782 79.510674 3.075769  76.857924
E0012.fasta 29  82  106 183 36  0.935934  64.199529 10.529236 6.406633  180.974063  71.992517 16.560459 2.44229 0.578834  0.054372  0.011758  0.035469  0.002909  59.380769 1.889886  0.11658 0.201405  0.908205  41.986265 77  24  126 0.435434  71.027563 12.470137 6.444319  181.581571  43.355201 7.05468 6.38779 180.701199  78.218863 45.45875  6.38779 180.701199
```
### Feature Importance (experimental)

Feature Importance using Random Forest-based non-linear regression between the features and GGI results

```Bash
ggpy fi -X demo/features_stats_demo.tsv -y demo/out_ggi_demo.txt
cat rf_FI_demo.csv
```
```
features        mean    var
invariants      0.023829008361325003    0.000492148638104501
gc_mean_pos3    0.013068959663030868    0.0002909344803794293
gap_var_pos3    0.005061315900347991    7.361446335854301e-05
nogap_prop      0.004502553450994113    8.56684923586634e-05
saturation      0.004119322498580811    0.0008387402170351275
seq_len_nogap   0.0036749280137925167   0.0002324035391180833
gap_mean        0.0035198531803458605   3.000429331750812e-05
gap_var_pos2    0.0034023410447559686   3.718083011942181e-05
gap_var_pos1    0.0033674991836586195   6.344403202970327e-05
gap_mean_pos1   0.003366633540142853    0.00011378203224143411
inter_len_mean  0.00273815317156398     0.0008308889860384133
inter_len_var   0.002459768056128803    0.005978171102088079
gap_mean_pos3   0.001839618363445573    1.3255972663730304e-05
gap_var 0.00072738945614525     2.7060906477610624e-05
gap_mean_pos2   0.00048259183547760674  1.2072408782394557e-05
pi_mean 0.0004722891538747615   0.00021575560888120713
patterns        0.000422681927520837    4.0679019633006905e-05
pi_std  8.654658290563468e-05   0.0005766083930162109
vars    7.418024850700645e-05   1.3976022285157172e-05
supp_mean       -0.0006585715613932451  0.0003747410242772111
total_tree_len  -0.0015306382925734353  0.0005201722323721086
seq_len -0.0018029590073972601  7.140911571362891e-05
gc_mean_pos2    -0.002237525762526324   0.002483625852550525
nheaders        -0.0034141967666995836  0.0006830433858946049
ter_len_mean    -0.0034281898822886205  0.00021689748391677419
singletons      -0.004800109089748769   0.00010866976775091261
treeness_o_rcv  -0.006688627586931427   0.00011069395128118889
gc_mean -0.006854598462143618   0.0007377179076071982
rcv     -0.006886856513351858   0.00013427768647595177
gc_var_pos3     -0.007278573565794007   9.523136628298002e-05
coeffVar_len    -0.007539838151231475   0.00016875566490389952
pis     -0.007605855262450881   6.907413350695042e-05
entropy -0.010740937363322168   7.580120884386808e-05
gc_var  -0.013119726647331848   0.00017455102520959055
ter_len_var     -0.017383212105414284   0.0009129408539352619
gc_var_pos1     -0.022562229308329795   0.0014302787795341092
gc_var_pos2     -0.025694528047872575   0.0008607445434747801
gc_mean_pos1    -0.02702702867996145    0.0023898373955544452
treeness        -0.029689025941821774   0.002322888932125225
LB_std  -0.09144113738811846    0.037928522325441
```



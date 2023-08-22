#!/usr/bin/env python

import argparse


# ITT_OUTFILENAME = "support_tt.tsv"
# RF_OUTFILENAME = "RF_distance.tsv"

parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                                    GGI and more
                                      ''')
subparsers = parser.add_subparsers(help='', dest='subcommand')

# ggi -------------------------------------------------------------------------
ggi = subparsers.add_parser('ggi',
                             help = "Gene-Genealogy Interrogation",
                             formatter_class = argparse.RawDescriptionHelpFormatter, 
                             description="""

                    Gene-Genealogy Interrogation

* Standard usage:

    $ %(prog)s [exon files] -t [taxonomy file]

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

    $ %(prog)s [exon files] -t [taxonomy file] -H [file with topologies]
        
        note: `-H` accepts a group-based topologies

* Specify pre-defined extended hypothesis:

    $ %(prog)s [exon files] -H [file with topologies] -e

        note: `-e` specifies that topologies at `-H` files are 
                extended (i.e., with actual species and not groups
                between parenthesis)

""")

ggi.add_argument('-H','--hypothesis',
                    metavar="",
                    default = None,
                    help='[Optional] Pre-defined hypothesis file, either extended or not')
ggi.add_argument('-e','--extended',
                    action="store_true",
                    help='''[Optional] If selected, file with topologies contains
                    extended trees (i.e., species instead of groups)''')

ggi.add_argument('-T','--parallel_gt',
                    action="store_true",
                    help='''[Optional] If selected, all cores are focused to one 
                    alignment at a time instead of one core per alignment 
                    (i.e., default strategy)''')

ggi.add_argument('-w','--write_extended',
                    action="store_true",
                    help='''[Optional] If selected, extended topologies are written and exit''')
ggi.add_argument('-r','--rooted',
                    action="store_true",
                    help='''[Optional] If selected, all posible rooted topologies
                        are generated when pre-defined hypothesis file is not given''')
ggi.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')
ggi.add_argument('-s','--suffix',
                    metavar="",
                    type= str,
                    default= "ggi.txt",
                    help='[Optional] Suffix for each written file [Default: ggi.txt]')


ggi_bl_raxml = ggi.add_argument_group('RAxML constrained tree parameters')
ggi_bl_raxml.add_argument('-E','--evomol',
                metavar="",
                type= str,
                default = 'GTRGAMMA',
                help='[Optional] RAxML evol. model for constrained tree inference [Default: GTRGAMMA]')
ggi_bl_raxml.add_argument('-c','--codon_aware',
                action="store_true",
                help='[Optional] If selected, codon partition file is added')
ggi_bl_raxml.add_argument('-i', '--iterations',
                metavar = "",
                type    = int,
                default = 1,
                help    = '[Optional] Number of iterations for MLEs [Default: 1]')


ggi_main_args = ggi.add_argument_group('required arguments')
ggi_main_args.add_argument('filenames',
                    metavar="alns",
                    nargs="*",
                    help='Alignments')   
ggi_main_args.add_argument('-t','--tax_file',
                    metavar="",
                    default = None,
                    # required= True,
                    help='Taxonomy file. Format in csv: "[sequence name],[group]"')    
ggi._action_groups.reverse()                     
# ggi -------------------------------------------------------------------------

# fstats -------------------------------------------------------------------------
fstats = subparsers.add_parser('features', help = "Features from both alignment and tree information",
                              formatter_class = argparse.RawDescriptionHelpFormatter, 
                              epilog = """
* codon_aware: currently only GC and Gap content
               can be obtained by codons

""",
                              description="""
                              
                    Summarize both alignment and tree information

Examples:

    * Standard usage:

        $ %(prog)s -A [alignment file extension] -T [tree file extension]

            Where '-A' and '-T' indicate file extensions for alignments and 
            trees, correspondingly. For each pair alignment and tree, extensions
            are taken out and leftovers should match

""")
fstats.add_argument('-A','--aln',
                    metavar="ext",
                    type= str,
                    # required=True,
                    default=".fasta",
                    help="Alignment file extension [Default = '.fasta']") 
fstats.add_argument('-T','--tree',
                    metavar="ext",
                    type= str,
                    # required=True,
                    default=".tree",
                    help="Tree file extension [Default = '.tree']")  
fstats.add_argument('-p','--path',
                    metavar="",
                    type= str,
                    default=".",
                    help="[Optional] Path to above files [Default: '.']")   
fstats.add_argument('-t','--taxonomy',
                    metavar="",
                    type= str,
                    required=False,
                    help='[Optional] Taxonomy file is CSV (e.g., "[seq],[Group]") [Default: None]')                      
# fstats.add_argument('-r','--reference',
#                      metavar="",
#                      type= str,
#                      default=None,
#                      help='[Optional] Reference tree [Default: None]')
# fstats.add_argument('-g','--group',
#                 metavar='',
#                 type= str,
#                 default=None,
#                 help='''[Optional] CSV-formated file containing alignmnet
#                         name and a label. This label is added as new column [Default: None]''')
fstats.add_argument('-c','--codon_aware',
                    action="store_true",
                    help='[Optional] If selected, metrics are obtained by codons*')
fstats.add_argument('-s','--suffix', 
                    metavar="",
                    type = str,
                    default='stats.tsv',
                    help='[Optional] suffix output names [Default = stats.tsv ]' )
fstats.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')          
# fstats -----------------------------------------------------------------------

# post_ggi --------------------------------------------------------------------------
post_ggi_sub = subparsers.add_parser('fi', help = "Feature importance using Random Forest",
                              formatter_class = argparse.RawDescriptionHelpFormatter, 
                              description="""

    Feature Importance (FI) using Random Forest-based non-linear regression 
               between the features and GGI results
            
Examples:

    * Standard usage:
        $ %(prog)s -X [feature matrix] -y [ggi results]
""")
post_ggi_sub.add_argument('-X','--features',
                      metavar = '',
                      required=True,
                      help='file from the "ggi" subcomand')
post_ggi_sub.add_argument('-y','--labels',
                      metavar = '',
                      required=True,
                      help='file from the "features" subcomand')
post_ggi_sub.add_argument('-t', '--test_size',
                      metavar = "",
                      type    = float,
                      default = 0.35,
                      help    = '[Optional] Testing size proportion [Default = 0.35]')
post_ggi_sub.add_argument('-e', '--trees',
                      metavar = "",
                      type    = int,
                      default = 100,
                      help    = '[Optional] number of estimators [Default = 100]')
post_ggi_sub.add_argument('-b', '--cv',
                      metavar = "",
                      type    = int,
                      default = 3,
                      help    = '[Optional] number of CV rounds for FI [Default = 3]')
post_ggi_sub.add_argument('-i', '--iter_perm',
                      metavar = "",
                      type    = int,
                      default = 30,
                      help    = '[Optional] number of iteration in permutation importance [Default = 30]')
post_ggi_sub.add_argument('-r', '--hyperparameter_size',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number hyperparameter tested in the specific range [Default = 1]') 
post_ggi_sub.add_argument('-s','--suffix', 
                        metavar="",
                        type = str,
                        default='demo',
                        help='[Optional] prefix name for outputs [Default = demo]' )
post_ggi_sub.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]') 
post_ggi_sub.add_argument('-p','--path',
                    metavar="",
                    type= str,
                    default=".",
                    help="[Optional] Path for outputs [Default: '.']")   
post_ggi_sub.add_argument('-g','--use_gini',
                    action="store_true",
                    help='[Optional] If selected, gini impurity metric is used to obtain feature importances')
# post_ggi --------------------------------------------------------------------------


def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "ggi":

        from ggpy.ggi import GGI

        GGI(
            sequences       = wholeargs.filenames,
            taxonomyfile    = wholeargs.tax_file,
            topologies      = wholeargs.hypothesis,
            are_extended    = wholeargs.extended,
            rooted          = wholeargs.rooted,
            codon_partition = wholeargs.codon_aware,
            threads         = wholeargs.threads,
            evomodel        = wholeargs.evomol,
            iterations      = wholeargs.iterations,
            parallel_gt     = wholeargs.parallel_gt,
            write_extended  = wholeargs.write_extended,
            suffix          = wholeargs.suffix,
        ).main()

    elif wholeargs.subcommand == "features":

        from fishlifetraits.stats import Features # expensive import

        Features(
            taxonomyfile   = wholeargs.taxonomy,
            path           = wholeargs.path,
            fasta_ext      = wholeargs.aln,
            tree_ext       = wholeargs.tree,
            reference_tree = None,
            groups_file    = None,
            codon_aware    = wholeargs.codon_aware,
            sym_tests      = False,
            threads        = wholeargs.threads,
            suffix         = wholeargs.suffix,
        ).write_stats()

    elif wholeargs.subcommand == "fi":
        # print(wholeargs)

        from ggpy.nonlinear_regression import FeatureImportance

        myclass = FeatureImportance(
            features_file    = wholeargs.features,
            all_ggi_results  = wholeargs.labels,
            n_jobs           = wholeargs.threads,
            n_estimators     = wholeargs.trees, #
            boots            = wholeargs.cv, #
            test_prop        = wholeargs.test_size, #
            n_repeats        = wholeargs.iter_perm, # iter for Permutations # 
            cv_size          = wholeargs.hyperparameter_size, # hyperparameter
            suffix           = wholeargs.suffix,
            out_folder       = wholeargs.path,
            gini_based       = wholeargs.use_gini # FI from gini impurity metric
        )
        myclass.regression_beta()
        # self.get_params()

if __name__ == "__main__":
    main()

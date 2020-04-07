import argparse
import subprocess
import numpy
from . import tcga_gather
from . import tcga_wrangle
from . import mcl_ea


def tcga_fi():
    # 6.a FI
    cm = "Rscript ../scripts/R/tcga_eda.R %s" % "tcga"
    subprocess.call(cm, shell=True)
    # sleep & check

    # 6.b mcl
    cm = "for i in `ls /opt/data/TCGA/...`; do ./cluster_mcl.sh $i tcga; done"
    subprocess.call(cm, shell=True)
    # sleep & check

    # 6.c EA
    mcl_ea.enrichment("tcga")
    # sleep & check


def tcga_adj():
    cm = "Rscript ../scripts/R/make_adj.R %s" % "tcga"
    subprocess.call(cm, shell=True)
    # sleep & check


def tcga_pipeline():
    """
    1. Gather TCGA sample files via The GDC API files end point https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/
    2. Use TCGAUtils package to download aliquotes ID's to deduplicate TCGA samples (ex. ran on multiple chips)
    3. Make project data frames csv files where rows are ENSG Ids and columns are samples (de-duplicated by step 2)
    4. Use BioMart to export unique-genes name to ENSG mapping via  https://uswest.ensembl.org/Help/Faq?id=125 most counts/mappings were fecthed via UI vs APIs
    5. EDA - QA/QC deduplicate genes, counts to cpm's, make gk_central & Functional Interactions (FIs) mappings in consortia
       tag and remove sample outlires via PCA
    6. FI correlations
        a. Make spearman correlation for FIs
        b. Weighted MCl clustering on FIs
        c. Reactome pathway enrichment analysis
    7. All genes correlations
        a. Make spearman correlation adjacency

    :return:
    """
    # 1. gather tcga
    tcga_gather.gather_samples()

    # 2. deduplicate TCGA samples
    cm = "Rscript ../scripts/R/tcga_clean_mapping.R"
    subprocess.call(cm, shell=True)
    # sleep & check

    # 3. make project data frames
    tcga_wrangle.build_matrix()
    # sleep & check

    # 4. BioMart uifetch mart_export-unique-genes.txt

    # 5. EDA - QA/QC
    cm = "Rscript ../scripts/R/tcga_eda.R"
    subprocess.call(cm, shell=True)
    # sleep & check

    # can sep each step in multiple processes
    # 6 FI
    # tcga_fi()

    # 7. all genes spearman Adj
    # tcga_adj()


def gtex_fi():
    # 3.a FI
    cm = "Rscript ../scripts/R/tcga_eda.R %s" % "gtex"
    subprocess.call(cm, shell=True)
    # sleep & check

    # 3.b mcl
    cm = "for i in `ls /opt/data/GTEx...`; do ./cluster_mcl.sh $i gtex; done"
    subprocess.call(cm, shell=True)
    # sleep & check

    # 3.c EA
    mcl_ea.enrichment("gtex")
    # sleep & check


def gtex_adj():
    cm = "Rscript ../scripts/R/make_adj.R %s" % "gtex"
    subprocess.call(cm, shell=True)
    # sleep & check


def gtex_pipeline():
    """
    1. Gather GTEx data wget script
    2. EDA - QA/QC deduplicate genes, counts to cpm's, make gk_central & Functional Interactions (FIs) mappings in consortia
       tag and remove sample outlires via PCA
    3. FI correlations
        a. Make spearman correlation for FIs
        b. Weighted MCl clustering on FIs
        c. Reactome pathway enrichment analysis
    4. All genes correlations
        a. Make spearman correlation adjacency

    :return:
    """
    # 1. gather gtex from https://www.gtexportal.org/home/datasets
    # cd dir
    # wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
    # unzip file

    # 2. EDA - QA/QC
    cm = "Rscript ../scripts/R/gtex_eda.R"
    subprocess.call(cm, shell=True)
    # sleep & check

    # 3. FI
    # gtex_fi()

    # 4. all genes spearman Adj
    # gtex_adj()


def buildParser():
    """
    Builds the argument parser for switch modules and evaluates results

    :return: switch modules output
    """
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(title='commands',
                                       description='The following commands are available:',
                                       help='For additional help: "switch <COMMAND> -h"')

    # switch commandline interface to objs
    parser_tcga = subparsers.add_parser('tcga', help='tcga pipeline')
    parser_tcga.set_defaults(func=tcga_pipeline)

    parser_tcga_fi = subparsers.add_parser('tcgaFI', help='tcga FI')
    parser_tcga_fi.set_defaults(func=tcga_fi)

    parser_tcga_adj = subparsers.add_parser('tcgaAdj', help='tcga Adj')
    parser_tcga_adj.set_defaults(func=tcga_adj)

    parser_gtex = subparsers.add_parser('gtex', help='gtex pipeline')
    parser_gtex.set_defaults(func=gtex_pipeline)

    parser_gtex_fi = subparsers.add_parser('gtexFI', help='gtex FI')
    parser_gtex_fi.set_defaults(func=gtex_fi)

    parser_gtex_adj = subparsers.add_parser('gtexAdj', help='gtex Adj')
    parser_gtex_adj.set_defaults(func=gtex_adj)

    return parser


def performMain(args):
    args.func(args)


def main():
    print("hello switch - let's start our data ops!")
    numpy.random.seed(1234)
    args = buildParser().parse_args()

    performMain(args)


if __name__ == "__main__":
    main()

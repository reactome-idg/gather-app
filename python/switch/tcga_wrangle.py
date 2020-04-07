import os
import numpy
import pandas


def _get_col(in_path, file, md, level):
    print(file)
    col = pandas.read_csv("".join([in_path, file]), sep='\t', header=None)
    col.columns = ['ensg', level]
    col = col[col.ensg.str.contains('ENSG')]
    sample_id = md.index[md.file_name.isin([file])][0]
    col.rename(columns={level: sample_id}, inplace=True)
    col.set_index(['ensg'], inplace=True)
    col.index = [s.split('.')[0] for s in col.index]

    return col


def build_matrix():
    """
    http://biit.cs.ut.ee/gprofiler/page/apis
    https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/307
    https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
    https://docs.gdc.cancer.gov/Encyclopedia/pages/HTSeq-Counts/
    https://docs.gdc.cancer.gov/Encyclopedia/pages/Harmonized_Data/
    https://docs.cancergenomicscloud.org/docs/tcga-grch38-metadata
    """
    level = 'count'
    out_path = "/opt/data/TCGA/processed/raw-count/"
    in_path = "/opt/data/TCGA/source/raw-data/"
    in_path_deduplicated_barcode = "/opt/data/TCGA/source/meta-data/count_deduplicated_barcode.csv"

    md = pandas.read_csv(in_path_deduplicated_barcode)
    md.set_index('submitter_id', inplace=True)
    projects = list(set(md['project.project_id']))

    ids = [list(md.aliquot_submitter_id[md['project.project_id'].isin([p])]) for p in projects]

    files = os.listdir(in_path)

    dat_list = []
    for p in projects:
        metadata = md[md['project.project_id'].isin([p])]
        project_files = [f for f in metadata['file_name'] if f in files]
        if project_files:
            col_list = [_get_col(in_path=in_path, file=f, md=md, level=level) for f in project_files]
            express_df = pandas.concat(col_list, axis=1, ignore_index=False)
            express_df.replace([numpy.inf, -numpy.inf, numpy.nan], 0)
            dat_list.append(express_df)
            express_df.to_csv("".join([out_path, p, "_COUNT.csv"]))

    pandas.DataFrame(express_df.index).to_csv("/opt/data/TCGA/source/ui-fetch/ensemble-fetch/ensg.csv", header=False, index=False)

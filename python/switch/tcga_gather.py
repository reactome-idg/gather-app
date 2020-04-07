import re
import os
import requests
from requests.exceptions import ConnectionError
import pandas
from pandas.io.json import json_normalize
import json
import gzip
import numpy
import random


def fetch():
    """
    https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/
    fetch metadata of interest - some metadata may be returned as a nested obj
    data model base
    https://docs.gdc.cancer.gov/Data_Submission_Portal/Users_Guide/Data_Submission_Walkthrough/

    df.shape # total cases 9760 out of 10642
    Note there are duplications in samples with multiple runs over multiple chips

    :return:
    """
    fields = [
        "access",
        "created_datetime",
        "data_category",
        "data_format",
        "data_type",
        "experimental_strategy",
        "file_id",
        "submitter_aliquot_ids",
        "file_name",
        "file_size",
        "file_state",
        "md5sum",
        "origin",
        "platform",
        "revision",
        "state",
        "submitter_id",
        "tags",
        "type",
        "analysis.input_files.platform",
        "analysis.workflow_type",
        "analysis.input_files.experimental_strategy",
        "analysis.metadata.read_groups.sequencing_center",
        "analysis.metadata.read_groups.sequencing_date",
        "analysis.metadata.read_groups.instrument_model",
        "analysis.metadata.read_groups.is_paired_end",
        "analysis.metadata.read_groups.library_name",
        "cases.case_id",
        "files.cases.aliquot_ids",
        "cases.submitter_aliquot_ids",
        "cases.submitter_analyte_ids",
        "cases.samples.sample_id",
        "cases.samples.sample_type",
        "cases.samples.is_ffpe",
        "cases.project.name",
        "cases.project.project_id",
        "cases.project.disease_type",
        "cases.project.primary_site",
        "cases.tissue_source_site.bcr_id",
        "cases.tissue_source_site.code",
        "cases.tissue_source_site.name",
        "cases.tissue_source_site.project",
        "cases.demographic.demographic_id",
        "cases.demographic.ethnicity",
        "cases.demographic.gender",
        "cases.demographic.race",
        "cases.demographic.state",
        "cases.demographic.submitter_id",
        "cases.demographic.updated_datetime",
        "cases.demographic.year_of_birth",
        "cases.demographic.year_of_death",
        "cases.diagnoses.age_at_diagnosis",
        "cases.diagnoses.classification_of_tumor",
        "cases.diagnoses.created_datetime",
        "cases.diagnoses.days_to_birth",
        "cases.diagnoses.days_to_death",
        "cases.diagnoses.days_to_last_follow_up",
        "cases.diagnoses.days_to_last_known_disease_status",
        "cases.diagnoses.days_to_recurrence",
        "cases.diagnoses.diagnosis_id",
        "cases.diagnoses.last_known_disease_status",
        "cases.diagnoses.morphology",
        "cases.diagnoses.primary_diagnosis",
        "cases.diagnoses.prior_malignancy",
        "cases.diagnoses.progression_or_recurrence",
        "cases.diagnoses.site_of_resection_or_biopsy",
        "cases.diagnoses.state",
        "cases.diagnoses.submitter_id",
        "cases.diagnoses.tissue_or_organ_of_origin",
        "cases.diagnoses.tumor_grade",
        "cases.diagnoses.tumor_stage",
        "cases.diagnoses.updated_datetime",
        "cases.diagnoses.vital_status",
    ]

    fields = ",".join(fields)

    files_endpt = "https://api.gdc.cancer.gov/files"

    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "files.cases.project.program.name",
                    "value": ["TCGA"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": ["RNA-Seq"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.analysis.workflow_type",
                    "value": ["HTSeq - Counts"]  # HTSeq - Counts, HTSeq - FPKM, HTSeq - FPKM-UQ is also available
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.cases.samples.sample_type",
                    "value": ["Primary Tumor"]  # only primary since additional primary means recurrent tumor (~mechanism)
                }
            }
        ]
    }

    params = {
        "filters": filters,
        "fields": fields,
        "format": "json",
        "size": "10000000"
    }

    headers = {
        "Content-Type": "application/json"
    }

    try:
        response = requests.post(files_endpt, headers=headers, json=params)
    except ConnectionError as e:
        print(e)

    if response.status_code == 200:
        r = response.json()

        json_dat = [j for j in r['data']['hits']]
        df = pandas.DataFrame(json_dat)

        return df
    else:
        print('In Fetch - status code returned a value of %s' % response.status_code)


def unnest_metadata():
    df = fetch()
    save = True
    metadata_filename = 'tcga-api-count-metadata.csv'

    analysis = df['analysis'].apply(pandas.Series)
    df['workflow_type'] = analysis['workflow_type']

    row_list = []
    for i, r in analysis.iterrows():
        row_list.append(json_normalize(r['input_files'][0]))
    input_files = pandas.concat(row_list, sort=False, ignore_index=True)
    df['platform'] = input_files['platform']

    row_list = []
    for i, r in df.iterrows():
        row_list.append(json_normalize(r['cases'][0]))
    cases = pandas.concat(row_list, sort=False, ignore_index=True)

    # filter by cases:
    # TCGA projects with annotated tissue and diagnoses metadata
    # Note not all samples have complete cases documented for diagnoses or identifiers
    # It is best to fetch the metadata from UI as the API calls don't fully render the search
    msk = cases['project.project_id'].str.contains('TCGA') & cases['project.primary_site'].notnull() & cases['diagnoses'].notnull()

    cases = cases[msk]
    df = df[msk]
    # unnest more ----------------------------------------------
    row_list = []
    for i, r in cases.iterrows():
        row_list.append(json_normalize(r['samples'][0]))
    samples = pandas.concat(row_list, sort=False, ignore_index=True)

    # incase null made it
    if cases['diagnoses'][0][0].keys():
        col_list = list(cases['diagnoses'][0][0].keys())

    row_list = []
    for i, r in cases.iterrows():
        if pandas.notnull(cases['diagnoses'][i]):
            row_list.append(json_normalize(r['diagnoses'][0]))
        else:
            row_list.append(pandas.DataFrame(columns=col_list))

    diagnoses = pandas.concat(row_list, sort=False, ignore_index=True)

    # bind info together  --------------------------------------
    cases = cases[['case_id', 'demographic.ethnicity', 'demographic.gender', 'demographic.race',
                   'demographic.state', 'demographic.submitter_id',
                   'demographic.updated_datetime', 'demographic.year_of_birth',
                   'demographic.year_of_death', 'project.disease_type',
                   'project.name', 'project.primary_site', 'project.project_id',
                   'tissue_source_site.bcr_id', 'tissue_source_site.code',
                   'tissue_source_site.name', 'tissue_source_site.project']]

    df = df[['file_id', 'file_name', 'file_size', 'id', 'md5sum', 'created_datetime',
             'data_category', 'data_format', 'data_type', 'experimental_strategy', 'submitter_id',
             'type', 'workflow_type', 'state', 'access']]

    sample_ids = [c.replace('_demographic', '') for c in cases['demographic.submitter_id']]

    cases.index = sample_ids
    df.index = sample_ids
    samples.index = sample_ids
    diagnoses.index = sample_ids

    final_df = pandas.concat([df, cases, samples, diagnoses], axis=1, ignore_index=False)
    final_df['uuid'] = [d[0].replace('_count', '') for d in final_df.submitter_id.iloc[:, 0:1].values]
    final_df.drop(columns=['submitter_id'], inplace=True)

    if save:
        final_df.to_csv("/opt/data/TCGA/source/meta-data/%s" % metadata_filename)

    return final_df


def _download(out_path, file_id):
    """
    download raw gzipped count data

    :param out_path:
    :param file_id:
    :return:
    """
    data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)

    try:
        response = requests.get(data_endpt, headers={"Content-Type": "application/json"})
    except ConnectionError as e:
        print(e)

    if response.status_code == 200:
        response_head_cd = response.headers["Content-Disposition"]
        file_name = re.findall("filename=(.+)", response_head_cd)[0]

        with open("".join([out_path, file_name]), "wb") as output_file:
            output_file.write(response.content)


def gather_samples():
    """
    For each file ID, fetched and unnested from TCGA-API, make's an API call
    :return: gzip file for each sample in "/opt/data/TCGA/source/raw-data/" directory
    """
    final_df = unnest_metadata()
    out_path = "/opt/data/TCGA/source/raw-data/"

    for i, r in final_df.iterrows():
        print(i, r['id'])
        _download(out_path, r['id'])

    # sometimes the api fetch halts prematurely due to network or tcga site issues
    # check if retry is needed for existing files
    file_names = set(final_df.file_name)
    files_downloaded = set(os.listdir(out_path))

    if len(file_names) != len(files_downloaded):
        f = file_names - files_downloaded
        f_ids = final_df.id[final_df.file_name.isin(f)]
        [_download(out_path, f_name) for f_name in f_ids]



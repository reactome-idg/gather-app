import pandas, numpy, os, reactome2py
from reactome2py import analysis, content, utils
from pandas.io.json import json_normalize
import scipy.stats as stats
import statsmodels.stats.multitest as multi


def _get_ea(f, in_path, out_path, low_path_df, min_cluster_size, all_genes):
    project = f.split('_')[0].replace('out.', '')
    print(project)

    df = pandas.read_csv("".join([in_path, f]), sep='\n', header=None, names=['gene_set'])
    df = df.apply(lambda x: x.str.replace('\t', ','))
    df['cluster_size'] = df.gene_set.map(lambda x: len(x.split(",")))
    df['gene_list'] = df.gene_set.map(lambda x: x.split(","))

    df = df[df.cluster_size > min_cluster_size]

    num_clusters = df.shape[0]
    print(num_clusters)
    num_pathway = low_path_df.shape[0]

    # create matrices for contingency table ------------------------------------------
    reactome_pathway_overlap = pandas.DataFrame({'all_genes': all_genes})

    for i, r in low_path_df.iterrows():
        reactome_pathway_overlap[r.stId] = reactome_pathway_overlap.all_genes.isin(r.genes)

    cluster_overlap = pandas.DataFrame({'all_genes': all_genes})

    for i, r in df.iterrows():
        cluster_overlap[df.index[i]] = cluster_overlap.all_genes.isin(r.gene_list)

    reactome_pathway_overlap.set_index('all_genes', inplace=True)
    cluster_overlap.set_index('all_genes', inplace=True)

    # -----------------------------------------------------------------------------
    # for each cluster, for each top/low pathway compute the enrichmnet analysis fishers exact test and
    # correct pathway p-values by all known genes count and compute rate of pathway significance over all known pathways
    # in this level -log10 adjpval we can also compute a rate by number of clusters

    ll = []

    for c, col1 in cluster_overlap.iteritems():
        l = []
        for p, col2 in reactome_pathway_overlap.iteritems():
            cross_table = pandas.crosstab(col1, col2)
            if cross_table.shape != (2, 2):
                cross_table = cross_table.append({True: 0, False: 0}, ignore_index=True)

            oddsratio, pvalue = stats.fisher_exact(cross_table, alternative='greater')
            overlap = "_".join([str(c), p])
            name = list(low_path_df.name[low_path_df.stId.isin([p])])[0]

            d = dict(overlap=overlap, cluster=str(c), stId=p, name=name, oddsratio=oddsratio, pvalue=pvalue)
            l.append(d)

        enrichment_df = pandas.DataFrame(l)
        corrected = multi.multipletests(enrichment_df.pvalue, alpha=0.05, method='fdr_bh', is_sorted=False,
                                        returnsorted=False)

        enrichment_df['fdr_bh'] = corrected[1]
        enrichment_df['sig_pathway'] = corrected[0]
        enrichment_df_sig = enrichment_df[enrichment_df.sig_pathway]

        enrichment_df['pathway_adjpvalue'] = -numpy.log10(enrichment_df_sig.fdr_bh) / num_pathway

        ll.append(enrichment_df)
    ea_df = pandas.concat(ll)

    df['cluster'] = df.index
    df.cluster = df.cluster.astype(int)
    ea_df.cluster = ea_df.cluster.astype(int)
    ea_df_all = pandas.merge(ea_df, df, how='left', on='cluster')

    ea_df_all.to_csv("".join([out_path, project, '_EA.csv']), index=False)


def enrichment(consortia):
    if consortia in ['TCGA', 'tcga']:
        in_path = "/opt/data/GTEx/processed/fi/mcl-spearman/"
        out_path = "/opt/data/GTEx/processed/fi/ea-spearman/"
    elif consortia in ['GTEx', 'GTEX', 'gtex']:
        in_path = "/opt/data/GTEx/processed/fi/mcl-spearman/"
        out_path = "/opt/data/GTEx/processed/fi/ea-spearman/"

    pathway_level_stIds = '/opt/data/curated/ELV_Pathway_Stable_IDs_Ver70.txt'
    elv = pandas.read_csv(pathway_level_stIds, header=None)

    gm = utils.gene_mappings()
    gm = json_normalize(gm)
    all_genes = list(set([a for b in gm.genes.values.tolist() for a in b]))

    low_path_df = gm[gm.stId.isin(list(elv[0].values))]
    min_cluster_size = 2

    files = os.listdir(in_path)
    [_get_ea(f, in_path, out_path, low_path_df, min_cluster_size, all_genes) for f in files]

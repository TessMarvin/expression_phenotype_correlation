#Author: Tess Marvin (tmarvin@nd.edu)
#Usage: python trait_phen.py
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from argparse import ArgumentParser
import numpy as np
import scipy.stats
import pingouin as pg

def sample_and_geneqc(expdata):
    '''
    Here, we'll remove poor samples and poor genes, and QC low read #s
    Low reads are any gene/sample reading < 5 counts; those are zeroed out
    Poor samples are defined as those with < 3000 genes with reads
    Poor genes are those that appear in < 20% of samples
    '''

    # Curate Samples
    #If there are less than 5 reads -- zero it out
    expdata = expdata.mask(expdata < 5, 0)
    #so np.count_nonzero with axis = 1 counts the number of nonzero counts for a gene across all columns (samples)
    #basically, this counts how many samples this gene is in
    genecounts = pd.Series(data = np.count_nonzero(expdata, axis = 1),
                             index = expdata.index)
    #so np.count_nonzero with axis = 0 counts the number of nonzero counts for a sample across all rows (genes)
    #basically this counts how many genes with reads each sample has
    samplecounts = pd.Series(np.count_nonzero(expdata,axis = 0),
                          index = expdata.columns)
    #this ensures that only the genes that are in more than 20% of samples are included
    goodgenes = genecounts[genecounts/samplecounts.size > 0.2]

    goodsamples = samplecounts[samplecounts > 3000]
    #the .loc function creates a dataframe that only contains the goodgenes and goodsamples
    allcur = expdata.loc[goodgenes.index, goodsamples.index]

    return allcur

def tpm_norm(expdata, probefile, fraglength = 250):
    '''
    TPM is (reads / effective gene length (in kb)) / (sample reads / 1e6)
    The elaborate nonsense below just calculates this, it's a bit of a
    mess but hey that's programming
    '''

    # Read in probe metadata, subtract frag length,
    # and convert to dict


    probeinfo = pd.read_csv(probefile, usecols = [1,6])
    probeinfo.index = probeinfo['Geneid']
    probeinfo.drop(columns = 'Geneid', inplace = True)
    probeinfo = probeinfo[probeinfo.index.isin(gctdata.index)]
    probeinfo['Length'] = probeinfo['Length'].apply(lambda x: x-fraglength)

    # Build output frame from input frame (cheating), calc length and
    # build lookup from samplename to total counts

    tpm = expdata.copy()
    curated_tpm = expdata.copy()

    # Iterate over rows- if non sample rows, just copy data, if sample
    # rows, calculate TPM

    # Note: it really doesn't like the way I did this and will throw
    # warnings- it does work, it's fine, don't worry about it, if you
    # are smarter than I am go ahead and rewrite it

    for (col, data) in expdata.iteritems():

        numreads = data.sum()

        tempframe = pd.concat([probeinfo, data], axis = 1)

        tempframe['rpkb'] = tempframe.iloc[:,1].divide(tempframe['Length'].divide(1000))
        tempframe['tpm'] = (tempframe['rpkb'] / data.sum()) * 1e6
        tempframe['tpm'] = tempframe['tpm'].clip(lower = 0)
        curated_tpm.loc[:,col] = tempframe['tpm']
    # Removing genes whose expression is ~0 in > 80% of samples
    goodgenes=[]
    for index, row in curated_tpm.iterrows():

        if (np.count_nonzero(row > 0.1) / len(row) > 0.2):
            goodgenes.append(index)

    curated_tpm = curated_tpm[curated_tpm.index.isin(goodgenes)]


    return curated_tpm

def fix_names(countdata):
    countdata = countdata.T
    '''
    This fixes the name specific issues in the NF54GFPxNHP4026 cross.
    The names have changed several times and been recorded in different
    formats, so I'll fix them with this.
    '''
    countdata.index = countdata.index.str.replace('\/','', regex = True)
    countdata.index = countdata.index.str.replace('ND5A5', 'AC075', regex = True)
    countdata.index = countdata.index.str.replace('ND6G8', 'AC125', regex = True)
    countdata.index = countdata.index.str.replace('N1', '', regex = True)
    countdata.index = countdata.index.str.replace('\\.', '', regex = True)
    countdata.index = countdata.index.str.replace('_4026', '_NHP4026', regex = True)
    countdata.index = countdata.index.str.replace('^4026', 'NHP4026', regex = True)
    countdata.index = countdata.index.str.replace('2H9', 'AC030', regex = True)
    countdata.index = countdata.index.str.replace('6E5', 'AC033', regex = True)
     # This is a hot mess of regular expression that just converts the entire
    # sample name to strain_##, where ## is the sampling timepoint. The formatting
    # is very inconsistent throughout, so a mess of replacements need to be carefully
    # made. This mess of substitutions makes them.
    countdata.index = countdata.index.str.replace('^GF_PL[\d]+[a,b,c]{0,1}_', '',regex = True)
    countdata.index = countdata.index.str.replace('[A,B,C]_', '', regex = True)
    countdata.index = countdata.index.str.replace('_[d]+$', '', regex = True)
    countdata.index = countdata.index.str.replace('_S.*', '', regex = True)
    countdata.index = countdata.index.str.replace('_[0-9]{3,4}$', '', regex = True)
    countdata.index = countdata.index.str.replace('hpi', '', regex = True)
    countdata.index = countdata.index.str.replace('_T', '_', regex = True)
    return countdata.T
def correlation_matrix(gene, df_correctnames, elo_data, cov4_data):
    '''

    '''
    #So I had to do this because some samples did not have any replicates
    #If that was the case list(df_correctnames.loc[gene, name_str]) resulted in error because list() cannot be used on single float
    how_many_reps = {}
    names4 = {}
    names30 = {}
    names44 = {}
    for name in df_correctnames.loc[gene].index:
        if name not in how_many_reps.keys():
            how_many_reps[name] = 1
        else:
            how_many_reps[name] += 1
    #For each progeny sample, save the TPM for the gene of interest
    #As most progeny have multiple replicates, save these as a list in a dictionary with progeny_name : list of TPM
    #The dictionary is specific to the hpi (4, 30, 44)
    for name in df_correctnames.loc[gene].index:
        name_str = str(name)
        #If the sample has _44 in its name, it is for the 44 hpi dictionary
        if '_44' in name_str:
        #If a previous index referenced this sample, skip it because we already have the list saved
            if name_str in names44.keys():
                continue
            #Else, save the list in the dictionary (save it in the way that works depending on if it is just one rep or more)
            else:
                if how_many_reps[name] == 1:
                    m = []
                    m.append(df_correctnames.loc[gene, name_str])
                    names44[name] = m
                else:
                    m = list(df_correctnames.loc[gene, name_str])
                    names44[name] = m
        #If the sample has _30 in its name, it is for the 30 hpi dictionary
        if '_30' in name:
            if name_str in names30.keys():
                continue
            else:
                if how_many_reps[name] == 1:
                    m = []
                    m.append(df_correctnames.loc[gene, name_str])
                    names30[name] = m
                else:
                    m = list(df_correctnames.loc[gene, name_str])
                    names30[name] = m
        #If the sample has _4 in its name, it is for the 4 hpi dictionary
        if '_4' in name:
            #We don't want to confuse 4 hpi and 44 hpi
            if '_44' in name:
                continue
            if name_str in names4.keys():
                continue
            else:
                if how_many_reps[name] == 1:
                    m = []
                    m.append(df_correctnames.loc[gene, name_str])
                    names4[name] = m
                else:
                    m = list(df_correctnames.loc[gene, name_str])
                    names4[name] = m
    x_vals = []
    y_vals = []
    dict_covs = {}
    progs = []
    for r in list(cov4_data.index):
        dict_covs[r] = []
    for progeny in names4.keys():
        name_of_prog = progeny.split("_")[0]
        progs.append(name_of_prog)
        if("NF54" in progeny):
            name_of_prog = name_of_prog.split("g")[0]
        if name_of_prog in elo_data.index and progeny in list(cov4_data):
            for i in list(names4[progeny]):
                x_vals.append(i)
                y_vals.append(float(elo_data.loc[name_of_prog]))
                for r in list(cov4_data.index):
                    if(progeny in list(cov4_data)):
                        dict_covs[r].append(cov4_data.loc[r, progeny])
                    else:
                        dict_covs[r].append(float("NaN"))

    r, p =scipy.stats.pearsonr(x_vals,y_vals)
    print(r)
    print(p)
    df = pd.DataFrame(list(zip(x_vals,y_vals)), columns = ['gene_expression', 'ELO'])
    df_covs = pd.DataFrame.from_dict(dict_covs)
    final_df = pd.concat([df,df_covs], axis=1)
    #['PL01', 'PL02', 'PL03', 'PL04', 'PL05', 'PL05a', 'PL05b', 'UNK', 'xcov', 'ycov', 'InferredCov1']
    corr_object = part_corr(data=final_df, x='gene_expression', y='ELO', covar = list(cov4_data.index), method = 'pearson')
    print(corr_object)
#THIS FUNCTION IS FROM PINGOUIN PACKAGE -- I had to modify it for my purposes -- I did not write this
#https://pingouin-stats.org/_modules/pingouin/correlation.html
#Author: Raphael Vallat
def part_corr(data=None, x=None, y=None, covar=None, x_covar=None,
                 y_covar=None, tail='two-sided', method='pearson'):
    from pingouin.utils import _flatten_list
    assert isinstance(data, pd.DataFrame), 'data must be a pandas DataFrame.'
    assert data.shape[0] > 2, 'Data must have at least 3 samples.'
    assert isinstance(x, (str, tuple)), 'x must be a string.'
    assert isinstance(y, (str, tuple)), 'y must be a string.'
    assert isinstance(covar, (str, list, type(None)))
    assert isinstance(x_covar, (str, list, type(None)))
    assert isinstance(y_covar, (str, list, type(None)))
    if covar is not None and (x_covar is not None or y_covar is not None):
        raise ValueError('Cannot specify both covar and {x,y}_covar.')
    assert x != covar, 'x and covar must be independent'
    assert y != covar, 'y and covar must be independent'
    assert x != y, 'x and y must be independent'
    # Check that columns exist
    col = _flatten_list([x, y, covar, x_covar, y_covar])
    if isinstance(covar, str):
        covar = [covar]
    if isinstance(x_covar, str):
        x_covar = [x_covar]
    if isinstance(y_covar, str):
        y_covar = [y_covar]

    assert all([c in data for c in col]), 'columns are not in dataframe.'
    # Check that columns are numeric
    assert all([data[c].dtype.kind in 'bfiu' for c in col])

    # Drop rows with NaN
    data = data[col].dropna()
    assert data.shape[0] > 2, 'Data must have at least 3 non-NAN samples.'

    # Standardize (= no need for an intercept in least-square regression) -- this does NOT work with my dummy variable for plates
    #So, only standardize for those variables that work
    for c in col:
        if(data[c].std(axis=0) != 0):
            data[c] = (data[c] - data[c].mean(axis=0)) / data[c].std(axis=0)
    if covar is not None:
        # PARTIAL CORRELATION
        cvar = np.atleast_2d(data[covar].to_numpy())
        beta_x = np.linalg.lstsq(cvar, data[x].to_numpy(), rcond=None)[0]
        beta_y = np.linalg.lstsq(cvar, data[y].to_numpy(), rcond=None)[0]
        res_x = data[x].to_numpy() - cvar @ beta_x
        res_y = data[y].to_numpy() - cvar @ beta_y
    else:
        # SEMI-PARTIAL CORRELATION
        # Initialize "fake" residuals
        res_x, res_y = data[x].to_numpy(), data[y].to_numpy()
        if x_covar is not None:
            cvar = np.atleast_2d(C[x_covar].to_numpy())
            beta_x = np.linalg.lstsq(cvar, C[x].to_numpy(), rcond=None)[0]
            res_x = C[x].to_numpy() - cvar @ beta_x
        if y_covar is not None:
            cvar = np.atleast_2d(C[y_covar].to_numpy())
            beta_y = np.linalg.lstsq(cvar, C[y].to_numpy(), rcond=None)[0]
            res_y = C[y].to_numpy() - cvar @ beta_y
    return pg.corr(res_x, res_y, method=method, tail=tail)
if __name__ == '__main__':
    probefile = './GENE2_probeinfo.csv'
    raw_couts_data= './GENE2_NOPROBE.count'
    elo_file = './ELO_NF54xART1.csv'
    t4_cov_file = './T4_qtl_1.PEER_covariates.txt'
    gctdata = pd.read_csv(raw_couts_data,
                          sep = '\t',
    					  index_col = 0)
    elo_data = pd.read_csv(elo_file, header = 0,index_col = 0)
    cov4_data= fix_names(pd.read_csv(t4_cov_file, sep = '\t', header = 0, index_col = 0))
    qc_counts_data= sample_and_geneqc(gctdata)
    normalized_counts= tpm_norm(qc_counts_data, probefile)
    normalized_fixed_names = fix_names(normalized_counts)
    for gene in normalized_fixed_names.index:
        correlation_matrix(gene,normalized_fixed_names, elo_data, cov4_data)

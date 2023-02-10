import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
from pybiomart import Server
import operator
import functools as ft
import random
from multiprocessing import Process, Queue, cpu_count
from scipy.stats import zscore
from scipy.stats import boxcox
from sklearn.preprocessing import MinMaxScaler
from itertools import compress
import math


pd.options.mode.chained_assignment = None  # default='warn'

R_SCRIPTS_DIR= '/PIPELINE/R_SCRIPTS/'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expr', help='expression file list', default=None, required=True)
    parser.add_argument('-i', '--idir', help='input dir', default=None, required=True)
    parser.add_argument('-o', '--odir', help='output dir', default=None, required=True)
    parser.add_argument('-m', '--meta', help='metadata file name', default=None, required=True)
    args = parser.parse_args()

    # initialize expression and meta data from files
    myrnaseqdata = RNASeqData(rnaSeqFileList=args.expr,
                              metaFileList=args.meta,
                              inputDir=args.idir,
                              outputDir=args.odir,
                              target_thresholds_per_metakey= None, # 'study',
                              target_threshold = 0.5,
                              sample_subsets_metakey = None, #'group',
                              sample_subsets= None, #['Flight', 'Ground', 'Vivarium'],
                              covariate_list= [], #['libprep'],
                              target='threshold',
                              env_key='env',
                              target_metakey_name= 'target', #'oro positivity (%)', #target
                              sample_key = 'subj_id', #'sample', #'subj_id'
                              feature_key = 'feature', #'gene', #'feature',
                              zero_count_percent = 0.8,
                              low_count_threshold = 50,
                              low_count_percent = 0.8,
                              high_mean_threshold = 0,
                              top_n_var = 0,
                              amplify = {'n':0, 'var':10, 'seed':23},
                              stack_xformations = True,
                              env='append', # env, append, xformation
                              verbose_R = False,
                              gene_filter= None, #'protein_coding',
                              filterFile = None, # '/Users/jcasalet/Desktop/RESEARCH/LIVER/DATA/JC/BIOMART/lipid-go-mart-export.tsv'
                              filterFileColumn = None, #'Gene_name'
                              xformations = ['merge_boxcox','merge_zscore','merge_std','merge_clr','merge_log','merge_sqrt'],
                              xform_all=None,
                              filter_mask=[], #['filterGenesByType']
                              norm_all=False,
                              filterCount=100000000,
                              splitInHalf = False,
                              callApplyFilter=False #True
                              )


    # vst can only be run on count data, NOT on zscores (and need to drop any meta cols from expr like env, threshold
    # use DESeq to shrink data with vsd, rld, and ntd transformations
    '''if not myrnaseqdata.vst is None:
        myrnaseqdata.shrinkExpression(myrnaseqdata.vst, myrnaseqdata.outputDir + '/expr_' + myrnaseqdata.vst + '.csv')

        myrnaseqdata.prep4Crisp(inputFile=myrnaseqdata.outputDir + '/expr_' + myrnaseqdata.vst + '.csv',
                                outputFileBase=myrnaseqdata.outputDir + '/' + myrnaseqdata.expr_outputfile_prefix + '_'  +
                                myrnaseqdata.vst)'''

    # prepare data for crisp consumption
    myrnaseqdata.prep4Crisp(inputFile=None,
                            outputFileBase=myrnaseqdata.outputDir + '/' + myrnaseqdata.expr_outputfile_prefix + '_crisp')


    # save meta for deseq2 now that threshold is calculated
    myrnaseqdata.save_meta(df=myrnaseqdata.metaDF, fileName=myrnaseqdata.outputDir + '/meta-4-deseq2.csv')

    # split in half for causalnex
    if myrnaseqdata.splitInHalf:
        df1, df2 = myrnaseqdata.splitStackInHalf(myrnaseqdata.xformation_stack['all'])
        myrnaseqdata.save_expr(inputDF=df1, fileName=myrnaseqdata.outputDir + '/expr_half_1.pkl')
        myrnaseqdata.save_expr(inputDF=df2, fileName=myrnaseqdata.outputDir + '/expr_half_2.pkl')


class RNASeqData():
    def __init__(self, inputDir ='.',
                 outputDir='.',
                 outputFilePrefix='expr',
                 rnaSeqFileList=None,
                 metaFileList=None,
                 target_thresholds_per_metakey=None,
                 target_threshold=0.5,
                 middle50_samples=False,
                 RScriptPath='/usr/local/bin/Rscript',
                 sample_subsets_metakey=None,
                 sample_subsets=None,
                 covariate_list=None,
                 target=None,
                 target_metakey_name=None,
                 sample_key = None,
                 env_key = None,
                 feature_key = None,
                 zero_count_percent=0.8,
                 low_count_threshold = 5,
                 low_count_percent = 0.8,
                 high_mean_threshold=0,
                 top_n_var = 0,
                 amplify={'n': 0, 'var': 10, 'seed': 23},
                 stack_xformations=False,
                 env=None,
                 verbose_R = False,
                 filterFile = None,
                 filterFileColumn = None,
                 gene_filter='protein_coding',
                 xformations = [],
                 xform_all = None,
                 filter_mask=None,
                 norm_all = True,
                 filterCount=0,
                 splitInHalf=False,
                 callApplyFilter=True
                 ):
        self.RScriptPath = RScriptPath
        self.inputDir = inputDir
        self.outputDir = outputDir
        self.expr_outputfile_prefix = outputFilePrefix
        self.target_thresholds_per_metakey=target_thresholds_per_metakey
        self.middle50_samples = middle50_samples
        self.covariate_list = covariate_list
        self.sample_subsets = sample_subsets
        self.sample_subsets_metakey = sample_subsets_metakey
        self.target = target
        self.target_metakey_name = target_metakey_name
        self.target_threshold = target_threshold
        self.sample_key = sample_key
        self.feature_key = feature_key
        self.env_key = env_key
        self.zero_count_percent = zero_count_percent
        self.low_count_threshold = low_count_threshold
        self.low_count_percent = low_count_percent
        self.high_mean_threshold = high_mean_threshold
        self.top_n_var = top_n_var
        self.amplify = amplify
        self.stack_xformations = stack_xformations
        self.xformation_stack = dict()
        self.meta_stack = dict()
        self.env = env
        self.verbose_R = verbose_R
        self.filterFile = filterFile
        self.filterFileColumn = filterFileColumn
        self.gene_filter = gene_filter
        self.xformations = xformations
        self.xform_all = xform_all
        self.filter_mask=filter_mask
        self.norm_all = norm_all
        self.filterCount = filterCount
        self.splitInHalf = splitInHalf
        self.callApplyFilter = callApplyFilter

        # set up metadata
        self.metaFileList = metaFileList.split(',')

        # set up expression df dict
        self.rnaSeqFileList = rnaSeqFileList.split(',')

        # first read all data into the expr and meta dicts and apply filters
        self.meta_dict = dict()
        self.rnaExprDataDict = dict()
        for e, m in zip(self.rnaSeqFileList, self.metaFileList):
            # read in expr data from file
            self.rnaExprDataDict[e] = pd.read_csv(self.inputDir + '/' + e, sep=',', header=0)
            # read in meta data from file
            self.meta_dict[e] = pd.read_csv(self.inputDir + '/' + m, sep=',', header=0)
            self.meta_dict[e].columns = list(map(lambda i: i.lower(), self.meta_dict[e].columns))
            # set all sample names in metadata to strings
            self.meta_dict[e][self.sample_key] = [str(x) for x in self.meta_dict[e][self.sample_key]]
            # set all sample names in expr data to strings
            self.rnaExprDataDict[e].columns = list([str(x) for x in self.rnaExprDataDict[e].columns])
            # filter out any samples from meta data that are not in sample_subsets
            if self.sample_subsets_metakey:
                self.meta_dict[e] = self.meta_dict[e][self.meta_dict[e][self.sample_subsets_metakey].isin(self.sample_subsets)]
            self.rnaExprDataDict[e] = self.transpose_df(self.rnaExprDataDict[e], cur_index_col=self.feature_key, new_index_col=self.sample_key)
            # filter out any samples from expr data that are not in meta_data
            self.rnaExprDataDict[e] = self.rnaExprDataDict[e][self.rnaExprDataDict[e][self.sample_key].isin(list(self.meta_dict[e][self.sample_key]))]
            # apply filters on each expr data in dict
            self.rnaExprDataDict[e], self.meta_dict[e] = self.apply_filter(self.rnaExprDataDict[e], self.meta_dict[e], mask=self.filter_mask)


        # test for NaN
        # df.columns[df.isna().any()].tolist()

        # combine all meta dicts into single df and save in file for deseq2
        self.metaDF = pd.concat(list(self.meta_dict.values()))

        # combine all expr dicts into single df and save in file for deseq2
        self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on=self.feature_key), list(self.rnaExprDataDict.values()))
        # permute samples
        self.expressionDF, self.metaDF = self.permuteSamples(expr_df=self.expressionDF, meta_df=self.metaDF, fileSuffix='_permuted')
        self.save_expr(inputDF=self.expressionDF,
                       fileName = self.outputDir + '/expr-4-deseq2.csv',
                       transpose=True,
                       cur_index_col=self.feature_key,
                       new_index_col=self.sample_key)
        # now save it for ica
        self.save_expr(inputDF=self.expressionDF, fileName = self.outputDir + '/expr-4-ica.csv')

        # set up samples

        if self.target_thresholds_per_metakey:
            self.target_thresholds_per_metakey_dict = dict()
            keep_samples = []
            for meta_val in set(self.metaDF[self.target_thresholds_per_metakey]):
                m50,f25,f50,g50,g75 = self.getMiddle50(meta_key=self.target_thresholds_per_metakey, meta_val=meta_val)
                print('f25 for ' + str(meta_val), str(f25))
                print('f50 for ' + str(meta_val), str(f50))
                print('g50 for ' + str(meta_val), str(g50))
                print('g75 for ' + str(meta_val), str(g75))

                if np.isnan(f25) and not np.isnan(g75):
                    self.target_thresholds_per_metakey_dict[meta_val] = g50
                elif np.isnan(g75) and not np.isnan(f25):
                    self.target_thresholds_per_metakey_dict[meta_val] = f50
                else:
                    self.target_thresholds_per_metakey_dict[meta_val] = (g50 + f50) / 2
                    # self.target_thresholds_per_metakey_dict[meta_val] = (g75 + f25) / 2
                print('threshold for ' + str(meta_val), str(self.target_thresholds_per_metakey_dict[meta_val]))
                keep_samples += m50


        # now add env and target to expressionDF
        if self.stack_xformations:
            self.build_stack_from_xformations()

    def build_stack_from_xformations(self):
        for xformation in self.xformations:
            eval(xformation)(obj=self)
        # combine all xformation_stack into individual expr matrices
        self.combine_xformation_stack()
        self.save_xformation_stack()
        self.save_meta_stack()
        self.xformation_stack['all'], self.meta_stack['all'] = self.permuteSamples(expr_df=self.xformation_stack['all'], meta_df=self.meta_stack['all'],
                            fileSuffix=str(self.xformations))

        # combine all xformation stacks into a single expr matrix
        self.xformation_stack['all'] = self.transpose_df(df=self.xformation_stack['all'], cur_index_col=self.feature_key,
                                                         new_index_col=self.sample_key)
        # and do some stdization over entire matrix
        if not self.xform_all is None:
            if self.xform_all is 'std':
                self.xformation_stack['all'] = self.my_stdizedf(df=self.xformation_stack['all'], across='features')
            elif self.xform_all is 'minmax':
                self.xformation_stack['all'] = self.my_minmaxscaler(df=self.xformation_stack['all'])
            elif self.xform_all is 'zscore':
                self.xformation_stack['all'] = self.my_zscore(df=self.xformation_stack['all'])
            elif self.xform_all is 'boxcox':
                self.xformation_stack['all'] = self.my_boxcox(df=self.xformation_stack['all'])
            elif self.xform_all is 'filterCount':
                samples = list(self.xformation_stack['all'][self.sample_key])
                df = self.xformation_stack['all'].drop(columns=[self.sample_key])
                df=df[df<self.filterCount]
                dropCols=pd.isnull(df).any(0)
                keepCols = [not i for i in dropCols]
                genes = list(compress(list(df.columns), keepCols))
                df[self.sample_key] = samples
                df=df[df.columns.intersection([self.sample_key] + genes)]
                self.xformation_stack['all'] = df
            else:
                print('unknown xformation for all: ', str(self.xform_all))



    def splitStackInHalf(self, df):
        thresh = list(df[self.target])
        env = list(df[self.env_key])
        df_t = self.transpose_df(df=df, cur_index_col=self.sample_key, new_index_col=self.feature_key)
        halfway_point = int(len(df_t)/2)
        df_1 = df_t[0:halfway_point]
        df_2 = df_t[halfway_point:]
        df_1_t = self.transpose_df(df=df_1, cur_index_col=self.feature_key, new_index_col=self.sample_key)
        df_2_t = self.transpose_df(df=df_2, cur_index_col=self.feature_key, new_index_col=self.sample_key)
        df_1_t[self.env_key] = env
        df_1_t[self.target] = thresh
        df_2_t[self.env_key] = env
        df_2_t[self.target] = thresh
        return df_1_t, df_2_t

    def save_meta_stack(self):
        for xformation in self.meta_stack:
            self.save_meta(self.meta_stack[xformation], fileName=self.outputDir + '/' + xformation + '-xformation-meta-stack.csv')

    def save_xformation_stack(self):
        for xformation in self.xformation_stack:
            self.save_expr(self.xformation_stack[xformation], fileName=self.outputDir + '/' + xformation + '-xformation-expr-stack.csv')

    def update_samples_in_xformation(self, xformation, raw=False):
        meta_df = pd.DataFrame()

        for i in range(len(self.xformation_stack[xformation])):
            old_name = self.xformation_stack[xformation].iloc[i][self.sample_key]
            new_name =  old_name + '__' + xformation
            self.xformation_stack[xformation].iloc[i, self.xformation_stack[xformation].columns.get_loc(self.sample_key)] = new_name
            meta_data = self.metaDF[self.metaDF[self.sample_key] == old_name]
            meta_data[self.sample_key] = new_name
            meta_data[self.env_key] = xformation
            meta_df = meta_df.append(meta_data)
        return self.xformation_stack[xformation], meta_df

    def combine_xformation_stack(self):

        # 1. tranpose to samples x genes
        # 2. change sample names per xformation
        # 3. create list of genes that intersect all xformations in stack
        counter=0
        for xformation in self.xformations:
            self.xformation_stack[xformation] = self.transpose_df(self.xformation_stack[xformation], cur_index_col=self.feature_key, new_index_col=self.sample_key)
            self.xformation_stack[xformation], self.meta_stack[xformation] = self.update_samples_in_xformation(xformation)
            if counter == 0:
                columns = self.xformation_stack[xformation].columns
                counter = 1
            else:
                columns = columns.intersection(self.xformation_stack[xformation].columns)
            print('len of cols for xformation ', xformation,  ' = ', str(len(self.xformation_stack[xformation].columns)))
            print('len of intersection = ', str(len(columns)))

        # 3. concatenate the expr from each xformation into one
        df = pd.concat([v for k,v in self.xformation_stack.items()])

        # 4. concatenate the meta from each xformation into one
        self.meta_stack['all'] = pd.concat([v for k,v in self.meta_stack.items()])

        # 5. intersect genes
        df = df[df.columns.intersection(list(columns))]

        # 6. transpose back to genes x samples
        self.xformation_stack['all'] = self.transpose_df(df, cur_index_col=self.sample_key, new_index_col=self.feature_key)


    def my_normalize(self, df, meta):
        inputFile = self.outputDir + '/expr_before_normalize.csv'
        outputFile = self.outputDir + '/expr_after_normalize.csv'
        metaFile = self.outputDir + '/temp-meta.csv'
        self.save_expr(df, inputFile)
        self.save_meta(meta, metaFile)
        cmd = ['/usr/local/bin/R', '-f', R_SCRIPTS_DIR + '/normalize.R', '--args', inputFile, metaFile, outputFile]
        self.callR(cmd)
        df = pd.read_csv(outputFile, sep=',', header=0)
        if self.feature_key in df.columns and 'Unnamed: 0' in df.columns:
            df.drop(columns=['Unnamed: 0'], inplace=True)
        return df

    def my_zscore(self, df):
        # this method requires df in samples x genes
        if self.sample_key in list(df.columns):
            transpose=False
            samples = list(df[self.sample_key])
            df1 = df.drop(columns=[self.sample_key])
        else:
            genes = list(df[self.feature_key])
            df1 = self.transpose_df(df, cur_index_col=self.sample_key, new_index_col=self.feature_key)
            df1 = df1.drop(columns=[self.feature_key])
            transpose=True
        df1 = df1.apply(zscore, axis=0)
        # self.expressionDF = (self.expressionDF - self.expressionDF.mean(axis=1)) / self.expressionDF.std(axis=1)
        if transpose:
            df1 = self.transpose_df(df1, cur_index_col=self.feature_key, new_index_col=self.sample_key)
            df1[self.feature_key] = genes
        else:
            df1[self.sample_key] = samples
        return df1

    def my_minmaxscaler(self, df):
        if self.feature_key in list(df.columns):
            df = self.transpose_df(df, cur_index_col=self.feature_key, new_index_col=self.sample_key)
            transpose = True
        else:
            transpose = False

        samples = list(df[self.sample_key])
        genes = list(df.columns)[1:]
        df1=df.drop(columns=[self.sample_key])
        scaler = MinMaxScaler()
        scaled_df = scaler.fit_transform(df1)
        scaled_df_pd = pd.DataFrame(scaled_df)
        scaled_df_pd.columns = genes
        scaled_df_pd[self.sample_key] = samples
        scaled_df_pd = scaled_df_pd[[self.sample_key] + genes]
        
        if transpose:
            return self.transpose_df(scaled_df_pd, cur_index_col=self.sample_key, new_index_col=self.feature_key)
        else:
            return scaled_df_pd

    def my_stdizedf(self, df, across='features'):
        if across == 'features':
            if self.feature_key in df.columns:
                features = list(df[self.feature_key])
                samples = list(df.columns)[1:]
                df1 = self.transpose_df(df, cur_index_col=self.feature_key, new_index_col=self.sample_key)
                df1 = df1.drop(columns=[self.sample_key])
                shape = 'feature_x_sample'
            elif self.sample_key in df.columns:
                samples = list(df[self.sample_key])
                features = list(df.columns)[1:]
                df1 = df.drop(columns=[self.sample_key])
                shape = 'sample_x_feature'

            df1 = (df1 - df1.mean()) / df1.std()
            df1[self.sample_key] = samples
            df1 = df1[[self.sample_key] + features]
            if shape == 'feature_x_sample':
                df1 = self.transpose_df(df1, cur_index_col=self.sample_key, new_index_col=self.feature_key)
            elif shape == 'sample_x_feature':
                pass
        elif across == 'samples':
            if self.feature_key in df.columns:
                features = list(df[self.feature_key])
                samples = list(df.columns)[1:]
                df1 = df.drop(columns=[self.feature_key])
                shape = 'feature_x_sample'
            elif self.sample_key in df.columns:
                samples = list(df[self.sample_key])
                features = list(df.columns)[1:]
                df1 = self.transpose_df(df, cur_index_col=self.sample_key, new_index_col=self.feature_key)
                df1 = df1.drop(columns=[self.feature_key])
                shape = 'sample_x_feature'

            df1 = (df1 - df1.mean()) / df1.std()
            df1[self.feature_key] = features
            df1 = df1[[self.feature_key] + samples]
            if shape == 'feature_x_sample':
                pass
            elif shape == 'sample_x_feature':
                df1 = self.transpose_df(df1, cur_index_col=self.feature_key, new_index_col=self.sample_key)

        return df1

    def my_sqrt(self, df=None):
        samples = list(df.columns)[1:]
        genes = list(df[self.feature_key])
        df1=df.drop(columns=[self.feature_key])
        sqrt_df = np.sqrt(df1)
        sqrt_df[self.feature_key] = genes
        sqrt_df = sqrt_df[[self.feature_key] + samples]
        sqrt_df = sqrt_df.reset_index().drop(columns=['index'])
        return sqrt_df

    def my_boxcox(self, df=None):
        # TODO parallelize me!
        print('dfmax = ', str(df.max()))
        if self.feature_key in list(df.columns):
            genes = list(df[self.feature_key])
            samples = list(df.columns[1:])
            df1=df.drop(columns=[self.feature_key])
            transpose = False
        else:
            samples = list(df[self.sample_key])
            genes = list(df.columns[1:])
            df1 = self.transpose_df(df, cur_index_col=self.sample_key, new_index_col=self.feature_key)
            df1 = df1.drop(columns=[self.feature_key])
            transpose = True

        df1 = df1 + 1
        bcDF = pd.DataFrame(columns=samples)
        for i in range(len(df1)):
            bc_array, bc_lambda = boxcox(df1.iloc[i])
            bcDF = bcDF.append(pd.DataFrame(np.array(bc_array).reshape(1, -1), columns=samples), ignore_index=True)

        bcDF[self.feature_key] = genes
        bcDF = bcDF[[self.feature_key] + samples]
        bcDF = bcDF.reset_index().drop(columns=['index'])

        print('bc lambda = ', str(bc_lambda))

        if transpose:
            return self.transpose_df(bcDF, cur_index_col=self.feature_key, new_index_col=self.sample_key)
        else:
            return bcDF


    def my_log_ratio(self, lr_type=None, df=None):
        import pyrolite
        # expr needs to be samples x genes and strictly positive
        samples = list(df.columns)[1:]
        genes = list(df[self.feature_key])
        df1=df.drop(columns=[self.feature_key])
        df1 = df1.T
        df1 = df1 + 1
        if lr_type == 'CLR':
            df1 = df1.pyrocomp.CLR()
        elif lr_type == 'ILR':
            df1 = df1.pyrocomp.ILR()
        df1 = df1.T
        df1[self.feature_key] = genes
        df1 = df1[[self.feature_key] + samples]
        df1 = df1.reset_index().drop(columns=['index'])
        return df1

    def my_log(self, base=10, expr=None):
        genes = list(expr[self.feature_key])
        df_np = np.array(expr.drop(columns=[self.feature_key]))
        if base == 10:
            df_np_log = pd.DataFrame(np.log10(df_np + 1))
        elif base == 2:
            df_np_log = pd.DataFrame(np.log2(df_np + 1))
        else:
            print('base unknown: ', base)
            sys.exit(1)
        df_np_log[self.feature_key] = genes
        cols = list(expr.columns[1:])
        df_np_log.columns = cols + [self.feature_key]
        df_np_log = df_np_log[[self.feature_key] + cols]
        return df_np_log

    def apply_filter(self, df, meta_df, mask=[]):
        # do all the filtering on raw data before normalizing or standardizing or log-xforming?
        # convert gene ids to gene names
        if self.callApplyFilter:
            if not 'convertIdsToNames' in mask:
                print('dims before converting to names: ', df.shape)
                df = self.convertIdsToNames(df=df)
                print('dims after converting to names: ', df.shape)

            # filter only protein-coding genes
            if not 'filterGenesByType' in mask:
                print('dims before filter by type: ', df.shape)
                if self.gene_filter:
                    df = self.filterGenesByType(df=df, gene_type=self.gene_filter, id=self.feature_key)
                print('dims after filter by type: ', df.shape)

            # filter only lipid-related genes (based on GO filters in biomart)
            if not 'filterGenesByFile' in mask:
                if self.filterFile is None or self.filterFileColumn is None:
                    pass
                else:
                    print('dims before filter by file: ', df.shape)
                    df = self.filterGenesByFile(df=df, fileName=self.filterFile, use=self.filterFileColumn)
                    print('dims after filter by file: ', df.shape)

            # filter genes with 0 count in at least p % samples
            if not 'filterGenesByPercentZeroCount' in mask:
                print('dims before filter 0: ', df.shape)
                df = self.filterGenesByPercentZeroCount(df, p=self.zero_count_percent)
                print('dims after filter 0: ', df.shape)

            # filter genes with count < n in at least p % samples
            if not 'filterGenesByPercentLowCount' in mask:
                print('dims before filter low: ', df.shape)
                df = self.filterGenesByPercentLowCount(df, n=self.low_count_threshold, p=self.low_count_percent)
                print('dims after filter low: ', df.shape)

            # filter genes with high means
            if not 'filterGenesByHighMeanThreshold' in mask:
                print('dims before filter high mean: ', df.shape)
                df = self.filterGenesByHighMeanThreshold(df, n=self.high_mean_threshold)
                print('dims after filter high mean: ', df.shape)
            # reduce number of genes to n top variance
            if not 'filterGenesByTopNSD' in mask:
                print('dims before filter by top n: ', df.shape)
                df= self.filterGenesByTopNSD(df, n=self.top_n_var)
                print('dims after filter by top n: ', df.shape)

            # amplify number of samples by n more samples
            if not 'amplify_expr' in mask:
                print('expr dims before amplify: ', df.shape)
                print('meta dims before amplify: ', meta_df.shape)
                df, meta_df = self.amplify_expr(df=df, metaDF=meta_df, amplify=self.amplify,  numProcs=4)
                print('dims after amplify: ', df.shape)
                print('meta dims after amplify: ', meta_df.shape)

        return self.transpose_df(df, cur_index_col=self.sample_key, new_index_col=self.feature_key), meta_df

    def setTargetByKey(self, key, exprDF, metaDF):
        valueSet = set(metaDF[key])
        if len(valueSet) != 2:
            print('need binary target')
            sys.exit(1)
        valueList = list(valueSet)
        targetDict = dict()
        for i in range(len(metaDF)):
            myKey = metaDF.iloc[i][self.sample_key]
            myValue = valueList.index(metaDF.iloc[i][key])
            targetDict[myKey] = myValue
        exprDF[key] = exprDF[self.sample_key].map(targetDict)

    def permuteSamples(self, expr_df, meta_df,fileSuffix=None):
        sample2index_dict = dict()
        print('meta sample keys: ', list(meta_df[self.sample_key]))
        print('expr sample keys: ', list(expr_df.columns))
        for i in range(meta_df.shape[0]):
            sample = str(meta_df.iloc[i][self.sample_key])
            j = expr_df.columns.get_loc(sample)
            sample2index_dict[sample] = (i, j)

        meta_df = meta_df.sample(frac=1)
        sample_list = list(meta_df[self.sample_key])
        expr_df = expr_df[[self.feature_key] + sample_list]
        self.expr_outputfile_prefix += '_permuted_'
        self.save_meta(meta_df, fileName=self.outputDir + '/meta-' + self.expr_outputfile_prefix + fileSuffix + '.csv')
        self.save_expr(inputDF=expr_df, fileName=self.outputDir + '/' + self.expr_outputfile_prefix + fileSuffix + '.csv')
        return expr_df, meta_df

    def prep4Crisp(self, inputFile, outputFileBase):
        if not inputFile is None:
            self.expressionDF = self.read_expr(inputFile)
        if self.stack_xformations:
            self.xformation_stack['all'] = self.add_env(covariate_list=self.covariate_list, expr_df=self.xformation_stack['all'], meta_df=self.meta_stack['all'], env=self.env)
        if self.target == 'threshold':
            if not inputFile is None:
                self.expressionDF = self.add_threshold(df=self.expressionDF, target_key = self.target_metakey_name,meta_df=self.metaDF, threshold=self.target_threshold)
            if self.stack_xformations:
                self.xformation_stack['all'] = self.add_threshold(df=self.xformation_stack['all'], target_key = self.target_metakey_name, meta_df=self.meta_stack['all'], threshold=self.target_threshold)
        else:
            #self.setTargetByKey(target)
            pass
        if not inputFile is None:
            self.save_expr(self.expressionDF, outputFileBase + '.csv')
            self.save_expr(self.expressionDF, outputFileBase + '.pkl')
        if self.stack_xformations:
            self.save_expr(self.xformation_stack['all'], outputFileBase + 'xformation_stack.csv')
            self.save_expr(self.xformation_stack['all'], outputFileBase + 'xformation_stack.pkl')

    def read_expr(self, fileName):
        df = pd.read_csv(fileName, header=0, sep=',', )
        first_col = df.columns[0]
        if first_col == self.feature_key:
            pass
        else:
            df.rename(columns = {first_col: self.feature_key}, inplace=True)

        return self.transpose_df(df, cur_index_col=self.feature_key, new_index_col=self.sample_key)

    def save_expr(self, inputDF=None, fileName=None, transpose=False, dropCols=[], cur_index_col=None, new_index_col=None):
        if inputDF is None:
            df = self.expressionDF
        else:
            df = inputDF
        if transpose:
            df = self.transpose_df(df.drop(columns=dropCols), cur_index_col, new_index_col)
            if fileName.endswith('.csv'):
                df.to_csv(fileName, sep=',', index=None)
            elif fileName.endswith('.pkl'):
                df.to_pickle(fileName)
            else:
                print('unknown filename extension: ', fileName)
                sys.exit(1)
        else:
            df = df.drop(columns=dropCols)
            if fileName.endswith('.csv'):
                df.to_csv(fileName, sep=',', index=None)
            elif fileName.endswith('.pkl'):
                df.to_pickle(fileName)
            else:
                print('unknown filename extension: ', fileName)
                sys.exit(1)

    def save_meta(self, df=None, fileName=None):
        if df is None:
            df = self.metaDF
        df.to_csv(fileName, sep=',', index=None)

    def convertIdsToNames(self, df):
        input_to_R = self.outputDir + '/expr_input.csv'
        output_from_R = self.outputDir + '/expr_output.csv'
        self.save_expr(inputDF=df, fileName = input_to_R, transpose=True, dropCols=[], cur_index_col=self.sample_key, new_index_col=self.feature_key)
        R_cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/convert_id_to_gene.R', '--args', input_to_R, output_from_R]
        self.callR(R_cmd)
        return self.read_expr(output_from_R)

    def shrinkExpression(self, transformation, fileName):
        if not transformation in ['vsd', 'rld', 'ntd']:
            print('transformation: ' + transformation + ' not known')
            sys.exit(1)
        input_to_R = self.outputDir + '/expr_before_shrink.csv'
        output_from_R = fileName
        meta_file = self.outputDir + '/meta.csv'
        self.save_expr(fileName = input_to_R, transpose=True, dropCols=[], cur_index_col=self.sample_key, new_index_col=self.feature_key)
        #self.save_expr(input_to_R)
        self.save_meta(fileName = meta_file)
        R_cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/shrink.R', '--args', input_to_R, meta_file, output_from_R,
                transformation, '--no-save']
        self.callR(R_cmd)

    def callR(self, cmd):
        import subprocess
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = proc.communicate(timeout=900)

        if self.verbose_R:
            print('Output: ' + o.decode('ascii'))
            print('Error: ' + str(e.decode('utf-8')))
            print('code: ' + str(proc.returncode))
        if str(proc.returncode) != '0':
            print('error in callR: exiting')
            sys.exit(1)

    def plotDict(self, ymin, ymax, myDict):
        fig, ax = plt.subplots()
        ax.set_ylim([ymin, ymax])
        ax.boxplot(myDict.values())
        ax.set_xticklabels(myDict.keys())
        fig.set_size_inches(12, 4)
        plt.show()

    def transpose_df(self, df, cur_index_col, new_index_col):
        if not cur_index_col in list(df.columns) or len(df) == 0:
            return df
        df = df.set_index(cur_index_col).T
        df.reset_index(level=0, inplace=True)
        cols = [new_index_col] + list(df.columns)[1:]
        df.columns = cols
        return df

    def add_threshold(self, df, meta_df, target_key, threshold=18):
        # create  threshold dictionary with arbitrary threshold
        threshold_dict = dict()
        threshold_dict_raw = dict()
        for i in range(len(meta_df)):
            key = meta_df.iloc[i][self.sample_key]
            key_raw = key.split('__')[0]
            val = meta_df.iloc[i][target_key]
            if self.target_thresholds_per_metakey:
                meta_key = meta_df.iloc[i][self.target_thresholds_per_metakey]
                thresh = self.target_thresholds_per_metakey_dict[meta_key]
            else:
                thresh = threshold
            if math.ceil(val) < thresh:
                threshold_dict[key] = 0
                threshold_dict_raw[key_raw] = 0
            else:
                threshold_dict[key] = 1
                threshold_dict_raw[key_raw] = 1

        # join val threshold dictionary to data frame
        df[self.target] = df[self.sample_key].map(threshold_dict)
        self.metaDF[self.target] = self.metaDF[self.sample_key].map(threshold_dict_raw)

        return df

    def create_env(self, covariate_list, meta_df):
        env_dict = dict()
        # create dictionary with sample id as key and concatenated column strings as value
        for i in range(len(meta_df)):
            key = meta_df.iloc[i][self.sample_key]
            counter = 0
            value = ''
            for e in covariate_list:
                if counter == 0:
                    value = str(meta_df.iloc[i][e])
                    counter += 1
                else:
                    value = value + ':' + str(meta_df.iloc[i][e])
            env_dict[key] = value
        return env_dict

    def add_env(self, covariate_list=['study', 'dissection', 'libprep'], expr_df=None, meta_df=None, env=None):
        env_dict = self.create_env(covariate_list, meta_df)

        # append or replace previously constructed env for xformation stack to env
        if not env is None:
            for i in range(len(expr_df)):
                sample = expr_df.iloc[i][self.sample_key]
                if self.env_key in meta_df.columns:
                    meta_env = meta_df[meta_df[self.sample_key] == sample][self.env_key].values[0]
                    if env == 'append':
                        env_dict[sample] += ':' + str(meta_env)
                    elif env == 'xformation':
                        env_dict[sample] = meta_env
                    elif env == self.env_key:
                        pass
                    elif env == 'xformation-append':
                        if 'clr' in meta_env or 'norm' in meta_env or 'std' in meta_env:
                            env_dict[sample] += ':' + str(meta_env)
                        else:
                            env_dict[sample] = ':' + str(env_dict[sample])

        # join env dictionary to data frame
        expr_df[self.env_key] = expr_df[self.sample_key].map(env_dict)
        return expr_df

    def filterGenesByHighMeanThreshold(self, df, n):
        if n == 0:
            pass
        else:
            genes = list(df.drop(columns=[self.sample_key]).columns)
            drop_genes = list()
            for gene in genes:
                if df[gene].mean() >= n:
                    drop_genes.append(gene)
            df = df.drop(columns=drop_genes)
        return df


    def filterGenesByPercentZeroCount(self, df=None, p=0):
        if p == 0:
            # this filtering assumes raw or normalized counts, not z-scores or log_xform
            pass
        else:
            df = self.transpose_df(df, self.sample_key, self.feature_key)
            df = df[(df == 0).sum(axis='columns') <= int(p * len(df.columns))]
            df = self.transpose_df(df, self.feature_key, self.sample_key)
        return df

    def filterGenesByPercentLowCount(self, df, n=0, p=0):
        if df is None:
            df = self.expressionDF
        if n == 0 or p == 0:
            # this filtering assumes raw or normalized counts, not z-scores or log_xform
            pass
        else:
            df = self.transpose_df(df, self.sample_key, self.feature_key)
            df = df[(df.loc[:, df.columns != self.feature_key] < n).sum(axis='columns') <= int(p * len(df.columns))]
            df= self.transpose_df(df, self.feature_key, self.sample_key)
        return df

    def filterGenesByFile(self, df, fileName, use=None):
        keepGenesDF = pd.read_csv(fileName, sep='\t', header=0)
        if not use in keepGenesDF.columns:
            print('use column name not in gene filter file: ', use)
            sys.exit(1)
        filter_genes = list(keepGenesDF[use])
        filter_columns = [self.sample_key] + filter_genes
        df = df[df.columns.intersection(filter_columns)]
        #df.columns = [self.sample_key] + filter_genes
        #df = df.loc[:, df.columns.notna()]
        return df

    def filterGenesByType(self, df, gene_type='protein_coding', id='id'):
        gene_types = {'ribozyme', 'protein_coding', 'rRNA', 'TEC', 'IG_D_pseudogene', 'snRNA', 'IG_LV_gene', 'pseudogene',
                      'IG_J_gene', 'transcribed_unitary_pseudogene', 'processed_pseudogene', 'IG_V_gene', 'Mt_tRNA',
                      'TR_J_pseudogene', 'miRNA', 'Mt_rRNA', 'sRNA', 'IG_C_pseudogene', 'IG_C_gene', 'TR_J_gene',
                      'IG_pseudogene', 'transcribed_processed_pseudogene', 'scRNA', 'lncRNA', 'TR_V_pseudogene',
                      'TR_V_gene', 'misc_RNA', 'TR_D_gene', 'translated_unprocessed_pseudogene',
                      'transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene',
                      'IG_V_pseudogene', 'scaRNA', 'TR_C_gene', 'IG_D_gene', 'snoRNA'}
        if not gene_type in gene_types:
            print('gene_type: ' + str(gene_type) + ' not recognized')
            sys.exit(1)
        server = Server(host='http://www.ensembl.org')
        dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['mmusculus_gene_ensembl'])
        gene_info = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype'])

        if id == 'id':
            filter_genes = list(gene_info[(gene_info['Gene type'] == gene_type) & (gene_info['Gene type'] != 'Mt_rRNA')]['Gene stable ID'])
            filter_columns = [self.sample_key] + filter_genes
            df = df[df.columns.intersection(filter_columns)]
            new_columns = list(df.drop(columns=[self.sample_key]))
            #gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene name'])
            gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene stable ID'])
            df.columns = [self.sample_key] + gene_names
        elif id == self.feature_key:
            filter_genes = list(gene_info[gene_info['Gene type'] == gene_type]['Gene name'])
            filter_columns = [self.sample_key] + filter_genes
            df = df[df.columns.intersection(filter_columns)]
        df = df.loc[:, df.columns.notna()]
        return df

    def filterGenesByTopNSD(self, df=None, n=0):
        # df is genes X samples
        # calculate var, sort cols into n highest vars, drop shape[1]-n cols
        # first find range of var and print to stdout
        if n == 0:
            pass
        else:
            df = self.transpose_df(df, self.sample_key, self.feature_key)
            #sdList = df.std(axis=1)
            sdList = df.var(axis=1)
            sdDict = {k: v for v, k in enumerate(sdList)}
            if n < 0:
                sdDictSorted = sorted(sdDict.items(), key=operator.itemgetter(0), reverse=False)
            else:
                sdDictSorted = sorted(sdDict.items(), key=operator.itemgetter(0), reverse=True)
            topN = sdDictSorted[0:abs(n)]
            indices = [x[1] for x in topN]
            df = df.iloc[indices]
            df = self.transpose_df(df, self.feature_key, self.sample_key)
        return df

    def collapseGeneCounts(self, exprDF):
        df = self.transpose_df(exprDF, self.sample_key, self.feature_key)
        dups = df[df.duplicated(self.feature_key, keep=False)].sort_values(self.feature_key)
        dups_cols = list(dups.columns)
        dups['index'] = list(dups.index)
        dups = dups[['index'] + dups_cols]
        indices = dict()
        for i in range(len(dups)):
            if dups.iloc[i][self.feature_key] in indices:
                indices[dups.iloc[i][self.feature_key]].append(dups.iloc[i]['index'])
            else:
                indices[dups.iloc[i][self.feature_key]] = [dups.iloc[i]['index']]
        collapsed_genes = dict()
        for gene in indices:
            collapsed_genes[gene] = df.iloc[indices[gene]].sum(axis=0, numeric_only=True)
            df.drop(indices[gene], axis=0,inplace=True)
        for gene in indices:
            df = df.append(collapsed_genes[gene], ignore_index=True)
            index=len(df)-1
            df.loc[index, self.feature_key] = gene
        exprDF = self.transpose_df(df, self.feature_key, self.sample_key)
        return exprDF


    def getMiddle50(self, meta_key, meta_val):

        f25 = float(self.metaDF[(self.metaDF[self.target_thresholds_per_metakey] == 'Flight') & (self.metaDF[meta_key] == meta_val)].describe().loc['25%'][self.target_metakey_name])
        f50 = float(self.metaDF[(self.metaDF[self.target_thresholds_per_metakey] == 'Flight') & (self.metaDF[meta_key] == meta_val)].describe().loc['50%'][self.target_metakey_name])
        g50 = float(self.metaDF[(self.metaDF[self.target_thresholds_per_metakey] != 'Flight') & (self.metaDF[meta_key] == meta_val)].describe().loc['50%'][self.target_metakey_name])
        g75 = float(self.metaDF[(self.metaDF[self.target_thresholds_per_metakey] != 'Flight') & (self.metaDF[meta_key] == meta_val)].describe().loc['75%'][self.target_metakey_name])

        group_df = self.metaDF[self.metaDF[meta_key] == meta_val]
        middle_group_samples = group_df[(group_df[self.target_metakey_name] >= g75) & (group_df[self.target_metakey_name] <= f25)][self.sample_key]

        return list(middle_group_samples), f25,f50,g50,g75

    def divide(self, n, d):
       res = list()
       qu = int(n/d)
       rm = n%d
       for i in range(d):
           if i < rm:
               res.append(qu + 1)
           else:
               res.append(qu)
       return res

    def getStartAndEnd(self, partitionSizes, threadID):
        start = 0
        for i in range(threadID):
            start += partitionSizes[i]
        end = start + partitionSizes[threadID]

        return start, end

    def process_samples_subset(self, q, samples, expr_samples_df, expr_df, meta_df, n, var, threadID, numProcs, seed, key):

        results = dict()
        partitionSizes = self.divide(len(samples), numProcs)
        start, end = self.getStartAndEnd(partitionSizes, threadID)
        #print('id: ', str(threadID), 'numProcs: ', str(numProcs), 'start: ', str(start), 'end: ', str(end))
        temp_meta_df = pd.DataFrame(columns=meta_df.columns)
        temp_expr_df = pd.DataFrame(columns=expr_df.columns)
        myRands = list()
        for i in range(n):
            for sample in samples[start:end]:
                # add new sample to expr data
                expr_row = expr_df[expr_df[self.sample_key] == sample].drop(columns=[self.sample_key])
                noise = np.random.normal(0, var, expr_row.shape)
                noised_expr_row = expr_row + noise
                noised_expr_row[noised_expr_row < 0 ] = 0
                myRand = random.randint(0, 1000000)
                while myRand in myRands:
                    myRand = random.randint(0, 1000000)
                myRands.append(myRand)
                new_sample = sample + '_' + str(threadID) + '_' + str(myRand)
                #self.samples.append(new_sample)
                #new_sample = sample + '_' + str((i+1) * n * threadID)

                noised_expr_row[self.sample_key] = new_sample
                temp_expr_df = temp_expr_df.append(noised_expr_row, ignore_index=False)

                # add new sample to meta data
                meta_row = meta_df[meta_df[key] == sample]
                new_meta_row = meta_row.copy(deep=True)
                new_meta_row[key] = new_sample
                temp_meta_df = temp_meta_df.append(new_meta_row, ignore_index=True)
        #return expr_df, meta_df
        results['expr_df'] = temp_expr_df
        results['meta_df'] = temp_meta_df
        q.put(results)

    def amplify_expr(self, df=None, metaDF = None, amplify = None,  numProcs=1):
        if amplify['n'] == 0:
            return df, metaDF

        expr_df = df
        meta_df = metaDF
        random.seed(amplify['seed'])

        # amplify set of samples that match column value
        q = Queue()
        processList = list()
        for i in range(numProcs):
            p = Process(target=self.process_samples_subset, args=(q,list(expr_df[self.sample_key]), expr_df[self.sample_key], expr_df,
                                                                  meta_df, amplify['n'], amplify['var'], i, numProcs, amplify['seed'], self.sample_key, ))
            p.start()
            processList.append(p)

        results = dict()
        for i in range(numProcs):
            results.update(q.get())
            expr_df = expr_df.append(results['expr_df'], ignore_index=False )
            meta_df = meta_df.append(results['meta_df'], ignore_index=False)

        for i in range(numProcs):
            print('joining thread: ', str(i))
            processList[i].join()

        self.expr_outputfile_prefix += '_amplify-' + str(amplify['n']) + '-' + str(amplify['var']) + '_'

        return expr_df, meta_df

def merge_std(obj):
    obj.xformation_stack['merge_std'] = obj.my_stdizedf(
        ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                  list(obj.rnaExprDataDict.values())), across='features')

def std_merge(obj):
    stdize_before_merge = dict()
    for e in obj.rnaExprDataDict:
        stdize_before_merge[e] = obj.my_stdizedf(obj.rnaExprDataDict[e], across='features')
    obj.xformation_stack['std_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                                                        list(stdize_before_merge.values()))

'''def stdsamples_merge(obj):
    stdize_before_merge = dict()
    for e in obj.rnaExprDataDict:
        stdize_before_merge[e] = obj.my_stdizedf(obj.rnaExprDataDict[e], across='samples')
    obj.xformation_stack['stdsamples_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=self.feature_key),
                                                  list(stdize_before_merge.values()))'''

def merge_stdsamples(obj):
    obj.xformation_stack['merge_stdsamples'] = obj.my_stdizedf(
        ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                  list(obj.rnaExprDataDict.values())), across='samples')

def zscore_merge(obj):
    zscore_before_merge = dict()
    for e in obj.rnaExprDataDict:
        zscore_before_merge[e] = obj.my_zscore(obj.rnaExprDataDict[e])
    obj.xformation_stack['zscore_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                                                  list(zscore_before_merge.values()))

def merge_zscore(obj):
    obj.xformation_stack['merge_zscore'] = obj.my_zscore(
        ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(obj.rnaExprDataDict.values())))

def merge_minmax(obj):
    obj.xformation_stack['merge_minmax'] = obj.my_minmaxscaler(
        ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(obj.rnaExprDataDict.values())))

def minmax_merge(obj):
    minmax_before_merge = dict()
    for e in obj.rnaExprDataDict:
        minmax_before_merge[e] = obj.my_minmaxscaler(obj.rnaExprDataDict[e])
    obj.xformation_stack['minmax_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                                                   list(minmax_before_merge.values()))

def merge_norm(obj):
    obj.xformation_stack['merge_norm'] = obj.my_normalize(
        ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(obj.rnaExprDataDict.values())),
        pd.concat(list(obj.meta_dict.values())))

def norm_merge(obj):
    norm_before_merge = dict()
    for e in obj.rnaExprDataDict:
        norm_before_merge[e] = obj.my_normalize(obj.rnaExprDataDict[e], obj.meta_dict[e])
    obj.xformation_stack['norm_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                                                           list(norm_before_merge.values()))

def merge_norm_std(obj):
    obj.xformation_stack['merge_norm_std'] = obj.my_stdizedf(
                                                obj.my_normalize(
                                                    ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                                                              list(obj.rnaExprDataDict.values())),
                                                    pd.concat(list(obj.meta_dict.values()))), across='features')

def norm_merge_std(obj):
    norm_before_merge_stdize = dict()
    for e in obj.rnaExprDataDict:
        norm_before_merge_stdize[e] = obj.my_normalize(obj.rnaExprDataDict[e], obj.meta_dict[e])
    obj.xformation_stack['norm_merge_std'] = obj.my_stdizedf(ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                                                                       list(norm_before_merge_stdize.values())), across='features')

'''def norm_std_merge(obj):
    norm_std_before_merge = dict()
    for e in obj.rnaExprDataDict:
        norm_std_before_merge[e] = obj.my_stdizedf(obj.my_normalize(obj.rnaExprDataDict[e], obj.meta_dict[e]), across='features')
    obj.xformation_stack['norm_std_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=self.feature_key),
                                                       list(norm_std_before_merge.values()))'''

def merge_log(obj):
    obj.xformation_stack['merge_log'] = obj.my_log(base=2, expr=ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(obj.rnaExprDataDict.values())))

def log_merge(obj):
    log_before_merge = dict()
    for e in obj.rnaExprDataDict:
        log_before_merge[e] = obj.my_log(base=2, expr=obj.rnaExprDataDict[e])
    obj.xformation_stack['log_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(log_before_merge.values()))

def merge_clr(obj):
    obj.xformation_stack['merge_clr'] = obj.my_log_ratio(lr_type='CLR', df=ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key),
                                                                                                 list(obj.rnaExprDataDict.values())))

def clr_merge(obj):
    clr_before_merge = dict()
    for e in obj.rnaExprDataDict:
        clr_before_merge[e] = obj.my_log_ratio('CLR', df=obj.rnaExprDataDict[e])
    obj.xformation_stack['clr_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(clr_before_merge.values()))

def merge_sqrt(obj):
    obj.xformation_stack['merge_sqrt'] = obj.my_sqrt(df=ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(obj.rnaExprDataDict.values())))

def sqrt_merge(obj):
    sqrt_before_merge = dict()
    for e in obj.rnaExprDataDict:
        sqrt_before_merge[e] = obj.my_sqrt(df=obj.rnaExprDataDict[e])
    obj.xformation_stack['sqrt_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(sqrt_before_merge.values()))

def merge_boxcox(obj):
    obj.xformation_stack['merge_boxcox'] = obj.my_boxcox(df=ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(obj.rnaExprDataDict.values())))

def boxcox_merge(obj):
    boxcox_before_merge = dict()
    for e in obj.rnaExprDataDict:
        boxcox_before_merge[e] = obj.my_sqrt(df=obj.my_boxcox(df=obj.rnaExprDataDict[e]))
    obj.xformation_stack['boxcox_merge'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(boxcox_before_merge.values()))

def raw(obj):
    obj.xformation_stack['raw'] = ft.reduce(lambda left, right: pd.merge(left, right, on=obj.feature_key), list(obj.rnaExprDataDict.values()))

if __name__ == "__main__":
    main()
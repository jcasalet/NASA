import subprocess
import pandas as pd
from statistics import mean, median
from sklearn import preprocessing
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
from pyrolite.util.synthetic import normal_frame, random_cov_matrix


pd.options.mode.chained_assignment = None  # default='warn'

R_SCRIPTS_DIR='/Users/jcasalet/Desktop/CODES/NASA/TRRAC/PIPELINE/R_SCRIPTS/'

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
                              normalize = 'afterMerge',
                              standardize = 'afterMerge',
                              log_xform=2,
                              vst=None,
                              oro_thresholds_per_study=True,
                              sample_subsets=['Flight', 'Ground', 'Vivarium', 'Basal'],
                              env_list= ['dataset', 'libprep'],
                              target='oro_thresh',
                              gene_filter = 'protein_coding',
                              zero_count_percent = 0.80,
                              low_count_threshold = 5,
                              low_count_percent = 0.80,
                              top_n_var = 20000,
                              amplify = {'n':0, 'var':10, 'seed':23},
                              stack_xformations = True
                              )

    # vst can only be run on count data, NOT on zscores (and need to drop any meta cols from expr like env, oro_thresh
    # use DESeq to shrink data with vsd, rld, and ntd transformations
    if not myrnaseqdata.vst is None:
        myrnaseqdata.shrinkExpression(myrnaseqdata.vst, myrnaseqdata.outputDir + '/expr_' + myrnaseqdata.vst + '.csv')
        myrnaseqdata.prep4Crisp(inputFile=myrnaseqdata.outputDir + '/expr_' + myrnaseqdata.vst + '.csv',
                                outputFileBase=myrnaseqdata.outputDir + '/' + myrnaseqdata.expr_outputfile_prefix + '_'  +
                                myrnaseqdata.vst,
                                env_list=myrnaseqdata.env_list,
                                target=myrnaseqdata.target)

    # prepare data for crisp consumption
    myrnaseqdata.prep4Crisp(inputFile=None,
                            outputFileBase=myrnaseqdata.outputDir + '/' + myrnaseqdata.expr_outputfile_prefix + '_crisp',
                            env_list=myrnaseqdata.env_list, target=myrnaseqdata.target)



class RNASeqData():
    def __init__(self, inputDir ='.',
                 outputDir='.',
                 outputFilePrefix='expr',
                 rnaSeqFileList=None,
                 metaFileList=None,
                 oro_scale='raw',
                 oro_thresholds_per_study=True,
                 middle50_samples=False,
                 RScriptPath='/usr/local/bin/Rscript',
                 normalize = None,
                 standardize=None,
                 vst=None,
                 log_xform=None,
                 sample_subsets=None,
                 env_list=None,
                 target=None,
                 gene_filter='protein_coding',
                 zero_count_percent=0.8,
                 low_count_threshold = 5,
                 low_count_percent = 0.8,
                 top_n_var = 0,
                 amplify={'n': 0, 'var': 10, 'seed': 23},
                 stack_xformations=False
                 ):
        self.RScriptPath = RScriptPath
        self.inputDir = inputDir
        self.outputDir = outputDir
        self.expr_outputfile_prefix = outputFilePrefix
        self.vst=vst
        self.oro_thresholds_per_study=oro_thresholds_per_study
        self.oro_scale=oro_scale
        self.middle50_samples = middle50_samples
        self.env_list = env_list
        self.normalize = normalize
        self.standardize = standardize
        self.log_xform = log_xform
        self.sample_subsets = sample_subsets
        self.env_list = env_list
        self.target = target
        self.gene_filter = gene_filter
        self.zero_count_percent = zero_count_percent
        self.low_count_threshold = low_count_threshold
        self.low_count_percent = low_count_percent
        self.top_n_var = top_n_var
        self.amplify = amplify
        self.stack_xformations = stack_xformations
        self.xformation_stack = dict()


        # append subsets to file name
        for group in self.sample_subsets:
            self.expr_outputfile_prefix += '-' + group
        self.expr_outputfile_prefix += '_'

        # set up metadata
        self.metaFileList = metaFileList.split(',')
        self.metaDict = dict()
        for f in self.metaFileList:
            self.metaDict[f] = pd.read_csv(self.inputDir + '/' + f, sep=',', header=0)
        self.metaDF = pd.concat(list(self.metaDict.values()), ignore_index=True)
        self.metaDF.columns = list(map(lambda i: i.lower(), self.metaDF.columns))
        print('meta dims before filtering out groups: ', self.metaDF.shape)
        self.metaDF = self.metaDF[self.metaDF['group'].isin(self.sample_subsets)]
        print('meta dims after filtering out groups: ', self.metaDF.shape)
        self.samples=list(self.metaDF['sample'])
        #self.save_meta(self.metaDF, self.outputDir + '/metadata.csv')

        # set up expression df dict
        self.rnaSeqFileList = rnaSeqFileList.split(',')



        # first read all data into the expr and meta dicts and apply filters
        self.meta_dict = dict()
        self.rnaExprDataDict = dict()
        for e, m in zip(self.rnaSeqFileList, self.metaFileList):
            self.rnaExprDataDict[e] = pd.read_csv(self.inputDir + '/' + e, sep=',', header=0)
            self.meta_dict[e] = pd.read_csv(self.inputDir + '/' + m, sep=',', header=0)
            self.rnaExprDataDict[e] = self.transpose_df(self.rnaExprDataDict[e], cur_index_col='gene', new_index_col='sample')
            self.rnaExprDataDict[e], self.meta_dict[e] = self.apply_filter(self.rnaExprDataDict[e], self.meta_dict[e], mask=[])

        # combine all meta dicts into single df
        self.metaDF = pd.concat(list(self.meta_dict.values()))
        self.metaDF.columns = list(map(lambda i: i.lower(), self.metaDF.columns))

        # case: normalize and optionally log-transform then standardize, all after merge
        # first, normalize
        if self.normalize == 'afterMerge':
            self.expr_outputfile_prefix += '_norm-after-merge_'
            if self.standardize == 'beforeMerge':
                # MOR doesn't work on negative count data
                print('incompatible: ', 'normalize:' + self.normalize, 'standardize: ' + self.standardize)
                sys.exit(1)
            # merge first
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            if self.stack_xformations:
                self.xformation_stack['raw-beforeMerge'] = self.expressionDF
            self.expressionDF = self.my_normalize(self.expressionDF, self.metaDF)
            if self.stack_xformations:
                self.xformation_stack['raw-beforeMerge_normalized-beforeMerge'] = self.expressionDF

            # second, log transform
            if not self.log_xform is None:
                if self.log_xform == 2 or self.log_xform == 10:
                    self.expressionDF = self.my_log(base=self.log_xform, expr=self.expressionDF)
                elif self.log_xform == 'ILR':
                    self.expressionDF = self.my_log_ratio(lr_type='ILR')
                elif self.log_xform == 'CLR':
                    self.expressionDF = self.my_log_ratio(lr_type='CLR')
                self.expr_outputfile_prefix += '_log-' + str(self.log_xform)
                if self.stack_xformations:
                    self.xformation_stack['raw-beforeMerge_normalized-beforeMerge_log-beforeMerge'] = self.expressionDF

            # third, standardize
            if self.standardize == 'afterMerge':
                self.expressionDF = self.my_standardize(self.expressionDF)
                if self.stack_xformations:
                    self.xformation_stack['raw-beforeMerge_normalized-beforeMerge_log-beforeMerge_stdize-afterMerge'] = self.expressionDF
                self.expr_outputfile_prefix += '_stdize-after-merge_'

        # case: normalize, log, standardize before merge
        elif self.normalize == 'beforeMerge':
            self.expr_outputfile_prefix += '_norm-before-merge_'
            self.meta_dict = dict()
            # first normalize
            for e, m in zip(self.rnaSeqFileList, self.metaFileList):
                self.rnaExprDataDict[e] = self.my_normalize(self.rnaExprDataDict[e], self.meta_dict[e])

            # second, log transforms
            if not self.log_xform is None:
                for e in self.rnaExprDataDict.keys():
                    self.rnaExprDataDict[e] = self.my_log(base=self.log_xform, expr=self.rnaExprDataDict[e])
                self.expr_outputfile_prefix += '_log' + str(self.log_xform) + '-before-merge_'

            # third, standardize
            if self.standardize == 'beforeMerge':
                for e in self.rnaExprDataDict.keys():
                    self.rnaExprDataDict[e] = self.my_standardize(self.rnaExprDataDict[e])
                self.expr_outputfile_prefix += '_stdize-before-merge_'

            # now merge all the data
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            if self.stack_xformations:
                self.xformation_stack['raw-beforeMerge'] = self.expressionDF
            # sub-case that norm before merge and then stdize after merge?
            if self.standardize == 'afterMerge':
                self.expressionDF = self.my_standardize(self.expressionDF)
                if self.stack_xformations:
                    self.xformation_stack['raw-beforeMerge_stdize-afterMerge']
                self.expr_outputfile_prefix += '_stdize-after-merge_'

            #self.metaDF = pd.concat(list(self.meta_dict.values()))
            #self.metaDF.columns = list(map(lambda i: i.lower(), self.metaDF.columns))


        # case: just standardize (no normalize) and potentially log before or after merge
        elif not self.standardize is None:

            if self.standardize == 'afterMerge':
                self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
                if self.stack_xformations:
                    self.xformation_stack['raw-afterMerge'] = self.expressionDF
                if not self.log_xform is None:
                    self.expressionDF = self.my_log(self.log_xform, self.expressionDF)
                    if self.stack_xformations:
                        self.xformation_stack['raw-afterMerge_fterMerge_log-afterMerge'] = self.expressionDF
                self.expressionDF = self.my_standardize(self.expressionDF)
                if self.stack_xformations:
                    self.xformation_stack['raw--afterMerge_stdize-afterMerge'] = self.expressionDF
                self.expr_outputfile_prefix += 'stdize-after-merge'

            elif self.standardize == 'beforeMerge':
                for e, m in zip(self.rnaSeqFileList, self.metaFileList):
                    if not self.log_xform is None:
                        self.rnaExprDataDict[e] = self.my_log(self.log_xform, self.rnaExprDataDict[e])
                    self.rnaExprDataDict[e] = self.my_standardize(self.rnaExprDataDict[e])
                self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
                if self.stack_xformations:
                    self.xformation_stack['stdize-beforeMerge'] = self.expressionDF
                #self.metaDF = pd.concat(list(self.meta_dict.values()))
                #self.metaDF.columns = list(map(lambda i: i.lower(), self.metaDF.columns))
                self.expr_outputfile_prefix += 'stdize-before-merge'

        # case: only log xform
        elif not self.log_xform is None:
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            if stack_xformations:
                self.xformation_stack['raw'] = self.expressionDF
            self.expressionDF = self.my_log(self.log_xform, self.expressionDF)
            if stack_xformations:
                self.xformation_stack['raw_log'] = self.expressionDF
            self.expr_outputfile_prefix += '_log' + str(self.log_xform) + '_'

        # case: no transformations, just merge and filter data
        elif self.normalize is None and self.standardize is None and self.log_xform is None:
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            if self.stack_xformations:
                self.xformation_stack['raw'] = self.expressionDF

        else:
            print('no such transformation combination value: ', 'normalize = ', self.normalize, 'standardize = ', self.standardize, 'log = ', self.log_xform)
            sys.exit(1)

        # combine all xformation_stack into a single expr matrix
        self.xformation_stack['all'] = self.combine_xformation_stack()

        self.save_xformation_stack()

        # create and save merged, untransformed data for deseq2??

        self.genes = list(self.expressionDF['gene'])

        # set up samples
        self.oro_thresholds_per_study_dict=dict()
        keep_samples = []
        if self.oro_thresholds_per_study:
            for group in set(self.metaDF['group']):
                for study in set(self.metaDF['study']):
                    m50,u25,l75 = self.getMiddle50('sample', group, study)
                    self.oro_thresholds_per_study_dict[study] = (u25 + l75) / 2
                    keep_samples += m50

        # subset expressionDF based on sample list
        if self.middle50_samples:
            self.expressionDF = self.expressionDF[self.expressionDF.columns.intersection(keep_samples)]
            #self.metaDF = self.metaDF[]

        # permute samples to be in same order in expr as in meta
        self.permuteSamples()

        # convert genes x samples to samples x genes
        self.expressionDF = self.transpose_df(df=self.expressionDF, cur_index_col='gene', new_index_col='sample')

        # combine same gene ids into one row
        print('dims before collapse: ', self.expressionDF.shape)
        self.collapseGeneCounts()
        print('dims after collapse: ', self.expressionDF.shape)

    def save_xformation_stack(self):
        for xformation in self.xformation_stack:
            self.save_expr(self.xformation_stack[xformation], fileName=self.outputDir + '/xformation-expr-stack.csv')

    def change_sampleNames_perXformation(self, xformation):
        for i in range(len(self.xformation_stack[xformation])):
            new_name = self.xformation_stack[xformation].iloc[i]['sample'] + '-' + xformation
            self.xformation_stack[xformation].iloc[i, self.xformation_stack[xformation].columns.get_loc('sample')] = new_name
        return self.xformation_stack[xformation]

    def combine_xformation_stack(self):
        # 1. tranpose to samples x genes
        # 2. change sample names per xformation
        # 3. create list of genes that intersect all xformations in stack
        counter=0
        for xformation in self.xformation_stack:
            self.xformation_stack[xformation] = self.transpose_df(self.xformation_stack[xformation], cur_index_col='gene', new_index_col='sample')
            self.xformation_stack[xformation] = self.change_sampleNames_perXformation(xformation)
            if counter == 0:
                columns = self.xformation_stack[xformation].columns
                counter = 1
            else:
                columns = columns.intersection(self.xformation_stack[xformation].columns)
            print('len of cols for xformation ', xformation,  ' = ', str(len(self.xformation_stack[xformation].columns)))
            print('len of intersection = ', str(len(columns)))

        # 3. concatenate the expr from each xformation into one
        df = pd.concat([v for k,v in self.xformation_stack.items()])

        # 4. intersect genes
        df = df[df.columns.intersection(list(columns))]

        # 5. transpose back to genes x samples
        return self.transpose_df(df, cur_index_col='sample', new_index_col='gene')


    def my_normalize(self, df, meta):
        inputFile = self.outputDir + '/expr_before_normalize.csv'
        outputFile = self.outputDir + '/expr_after_normalize.csv'
        metaFile = self.outputDir + '/temp-meta.csv'
        self.save_expr(df, inputFile)
        self.save_meta(meta, metaFile)
        cmd = ['/usr/local/bin/R', '-f', R_SCRIPTS_DIR + '/normalize.R', '--args', inputFile, metaFile, outputFile]
        self.callR(cmd)
        self.expressionDF = pd.read_csv(outputFile, sep=',', header=0)
        if 'gene' in df.columns and 'Unnamed: 0' in df.columns:
            df.drop(columns=['Unnamed: 0'], inplace=True)
        return df

    def my_standardize(self, df):
        # this method requires df in samples x genes
        if 'gene' in list(df.columns):
            transpose=True
            df = self.transpose_df(df, cur_index_col='gene', new_index_col='sample')
        numeric_cols = list(df.columns[1:])
        samples = list(df['sample'])
        df.drop(columns=['sample'], inplace=True)
        df = df.apply(zscore, axis=0)
        # self.expressionDF = (self.expressionDF - self.expressionDF.mean(axis=1)) / self.expressionDF.std(axis=1)
        df['sample'] = samples
        df = df[['sample'] + numeric_cols]
        if transpose:
            return self.transpose_df(df, cur_index_col='sample', new_index_col='gene')
        else:
            return df

    def apply_filter(self, df, meta_df, mask=[]):
        # do all the filtering on raw data before normalizing or standardizing or log-xforming?
        # convert gene ids to gene names
        if not 'convertIdsToNames' in mask:
            print('dims before converting to names: ', df.shape)
            df = self.convertIdsToNames(df=df)
            print('dims after converting to names: ', df.shape)

        # filter only protein-coding genes
        if not 'filterGenesByType' in mask:
            print('dims before filter by type: ', df.shape)
            df = self.filterGenesByType(df=df, gene_type=self.gene_filter, id='gene')
            print('dims after filter by type: ', df.shape)

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

        # reduce number of genes to n top variance
        if not 'filterGenesByTopNSD' in mask:
            print('dims before filter by top n: ', df.shape)
            df= self.filterGenesByTopNSD(df, n=self.top_n_var)
            print('dims after filter by top n: ', df.shape)

        # amplify number of samples by n more samples
        if not 'amplify_expr' in mask:
            print('expr dims before amplify: ', df.shape)
            print('meta dims before amplify: ', meta_df.shape)
            df, meta_df = self.amplify_expr(df=df, metaDF=meta_df, amplify=self.amplify, key='sample', numProcs=4)
            print('dims after amplify: ', df.shape)
            print('meta dims after amplify: ', meta_df.shape)


        return self.transpose_df(df, cur_index_col='sample', new_index_col='gene'), meta_df

    def my_log_ratio(self, lr_type=None):
        import pyrolite
        # expr needs to be samples x genes and strictly positive
        samples = list(self.expressionDF.columns)[1:]
        genes = list(self.expressionDF['gene'])
        df=self.expressionDF.drop(columns=['gene'])
        df=df.T
        df=df+1
        if lr_type=='CLR':
            df = df.pyrocomp.CLR()
        elif lr_type == 'ILR':
            df = df.pyrocomp.ILR()
        df = df.T
        df['gene'] = genes
        df = df[['gene'] + samples]
        df = df.reset_index().drop(columns=['index'])
        return df

    def my_log(self, base=10, expr=None):
        genes = list(expr['gene'])
        df_np = np.array(expr.drop(columns=['gene']))
        if base == 10:
            df_np_log = pd.DataFrame(np.log10(df_np+1))
        elif base == 2:
            df_np_log = pd.DataFrame(np.log2(df_np+1))
        else:
            print('base unknown: ', base)
            sys.exit(1)
        df_np_log['gene'] = genes
        cols = list(expr.columns[1:])
        df_np_log.columns = cols + ['gene']
        df_np_log = df_np_log[['gene'] + cols]
        return df_np_log


    def setTargetByKey(self, key):
        valueSet = set(self.metaDF[key])
        if len(valueSet) != 2:
            print('need binary target')
            sys.exit(1)
        valueList = list(valueSet)
        targetDict = dict()
        for i in range(len(self.metaDF)):
            myKey = self.metaDF.iloc[i]['sample']
            myValue = valueList.index(self.metaDF.iloc[i][key])
            targetDict[myKey] = myValue
        self.expressionDF[key] = self.expressionDF['sample'].map(targetDict)

    def permuteSamples(self, metaKey='sample', exprKey='gene'):
        sample2index_dict = dict()
        df = self.expressionDF
        for i in range(self.metaDF.shape[0]):
            sample = self.metaDF.iloc[i][metaKey]
            j = df.columns.get_loc(sample)
            sample2index_dict[sample] = (i, j)

        self.metaDF = self.metaDF.sample(frac=1)
        sample_list = list(self.metaDF[metaKey])
        self.expressionDF = df[[exprKey] + sample_list]
        self.expr_outputfile_prefix += '_permuted_'
        self.save_meta(self.metaDF, fileName=self.outputDir + '/meta-' + self.expr_outputfile_prefix + '.csv')
        self.save_expr(inputDF=self.expressionDF, fileName=self.outputDir + '/' + self.expr_outputfile_prefix + '.csv')

    def prep4Crisp(self, inputFile, outputFileBase, env_list, target):
        if not inputFile is None:
            self.expressionDF = self.read_expr(inputFile)
        self.expressionDF = self.add_env(env_list=env_list, expr_df=self.expressionDF, meta_df=self.meta_Df)
        #self.xformation_stack['all'] = self.add_env(env_list=env_list, expr_df=self.xformation_stack['all'], meta_df=self.meta_stack)
        if target == 'oro_thresh':
            self.add_oro()
        else:
            self.setTargetByKey(target)
        self.save_expr(self.expressionDF, outputFileBase + '.csv')
        self.save_expr(self.expressionDF, outputFileBase + '.pkl')

    def mergeExprData(self, dataDict):
        self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(dataDict.values()))

    def read_expr(self, fileName):
        df = pd.read_csv(fileName, header=0, sep=',', )
        first_col = df.columns[0]
        if first_col == 'gene':
            pass
        else:
            df.rename(columns = {first_col: 'gene'}, inplace=True)

        return self.transpose_df(df, cur_index_col='gene', new_index_col='sample')

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
        self.save_expr(inputDF=df, fileName = input_to_R, transpose=True, dropCols=[], cur_index_col='sample', new_index_col='gene')
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
        self.save_expr(fileName = input_to_R, transpose=True, dropCols=[], cur_index_col='sample', new_index_col='gene')
        #self.save_expr(input_to_R)
        self.save_meta(fileName = meta_file)
        R_cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/shrink.R', '--args', input_to_R, meta_file, output_from_R,
                transformation, '--no-save']
        self.callR(R_cmd)

    def callR(self, cmd):
        import subprocess
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        o, e = proc.communicate(timeout=900)

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

    def add_oro(self, threshold=18):
        # create oro threshold dictionary with arbitrary threshold
        # >18 = 0 which is the 70% oro value, and there are 70% gc and 30% flight
        # or use 50% (median) oro value to be un-biased (13.96)
        # threshold=13.96 # median oro value
        # threshold=18 #70% oro value
        orothresh_dict = dict()
        for i in range(len(self.metaDF)):
            key = self.metaDF.iloc[i]['sample']
            oro = self.metaDF.iloc[i]['oro positivity (%)']
            study = self.metaDF.iloc[i]['study']
            if not self.oro_thresholds_per_study:
                thresh = self.oro_thresholds_per_study_dict[study]
            else:
                thresh = threshold
            if oro < thresh:
                orothresh_dict[key] = 0
            else:
                orothresh_dict[key] = 1

        # join oro threshold dictionary to data frame
        self.expressionDF['oro_thresh'] = self.expressionDF['sample'].map(orothresh_dict)

    def transpose_df(self, df, cur_index_col, new_index_col):
        df = df.set_index(cur_index_col).T
        df.reset_index(level=0, inplace=True)
        cols = [new_index_col] + list(df.columns)[1:]
        df.columns=cols
        return df

    def add_env(self, env_list=None, expr_df=None, meta_df=None):
        # create dictionary with sample id as key and concatenated column strings as value
        env_dict = dict()
        for i in range(len(meta_df)):
            key = meta_df.iloc[i]['sample']
            if not env_list:
                value = meta_df.iloc[i]['study'] + ':' + meta_df.iloc[i]['dissection'] + ':' + meta_df.iloc[i]['libprep']
            else:
                counter = 0
                for e in env_list:
                    if counter == 0:
                        value = str(meta_df.iloc[i][e])
                        counter+=1
                    else:
                        value = value + ':' + str(meta_df.iloc[i][e])
            env_dict[key] = value

        # join env dictionary to data frame
        expr_df['env'] = expr_df['sample'].map(env_dict)

        return expr_df

    def filterGenesByPercentZeroCount(self, df=None, p=0):
        if p == 0:
            # this filtering assumes raw or normalized counts, not z-scores or log_xform
            pass
        else:
            df = self.transpose_df(df, 'sample', 'gene')
            df = df[(df == 0).sum(axis='columns') <= int(p * len(df.columns))]
            df = self.transpose_df(df, 'gene', 'sample')
        return df

    def filterGenesByPercentLowCount(self, df, n=0, p=0):
        if df is None:
            df = self.expressionDF
        if n == 0 or p == 0:
            # this filtering assumes raw or normalized counts, not z-scores or log_xform
            pass
        else:
            df = self.transpose_df(df, 'sample', 'gene')
            df = df[(df.loc[:, df.columns != 'gene'] < n).sum(axis='columns') <= int(p * len(df.columns))]
            df= self.transpose_df(df, 'gene', 'sample')
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
            filter_genes = list(gene_info[gene_info['Gene type'] == gene_type]['Gene stable ID'])
            filter_columns = ['sample'] + filter_genes
            df = df[df.columns.intersection(filter_columns)]
            new_columns = list(df.drop(columns=['sample']))
            #gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene name'])
            gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene stable ID'])
            df.columns = ['sample'] + gene_names
        elif id == 'gene':
            filter_genes = list(gene_info[gene_info['Gene type'] == gene_type]['Gene name'])
            filter_columns = ['sample'] + filter_genes
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
            df = self.transpose_df(df, 'sample', 'gene')
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
            df = self.transpose_df(df, 'gene', 'sample')
        return df

    def collapseGeneCounts(self):
        df = self.transpose_df(self.expressionDF, 'sample', 'gene')
        dups = df[df.duplicated('gene', keep=False)].sort_values('gene')
        dups_cols = list(dups.columns)
        dups['index'] = list(dups.index)
        dups = dups[['index'] + dups_cols]
        indices = dict()
        for i in range(len(dups)):
            if dups.iloc[i]['gene'] in indices:
                indices[dups.iloc[i]['gene']].append(dups.iloc[i]['index'])
            else:
                indices[dups.iloc[i]['gene']] = [dups.iloc[i]['index']]
        collapsed_genes = dict()
        for gene in indices:
            collapsed_genes[gene] = df.iloc[indices[gene]].sum(axis=0, numeric_only=True)
            df.drop(indices[gene], axis=0,inplace=True)
        for gene in indices:
            df = df.append(collapsed_genes[gene], ignore_index=True)
            index=len(df)-1
            df.loc[index, 'gene'] = gene
        self.expressionDF = self.transpose_df(df, 'gene', 'sample')


    def getMiddle50(self, key, group, study):
        # group_25=float(df[(df['Group']==group) & (df['Study']==study)].describe().loc['25%']['ORO Positivity (%)'])
        # group_75=float(df[(df['Group']==group) & (df['Study']==study)].describe().loc['75%']['ORO Positivity (%)'])
        if len(self.metaDF[(self.metaDF['group'] == group) & (self.metaDF['study'] == study)]) == 0:
            return [], 0, 0
        group_25 = np.quantile(self.metaDF[(self.metaDF['group'] == group) & (self.metaDF['study'] == study)]['oro positivity (%)'], q=0.25)
        group_75 = np.quantile(self.metaDF[(self.metaDF['group'] == group) & (self.metaDF['study'] == study)]['oro positivity (%)'], q=0.75)
        group_df = self.metaDF[(self.metaDF['group'] == group) & (self.metaDF['study'] == study)]

        middle_group_samples = group_df[(group_df['oro positivity (%)'] >= group_25) & (group_df['oro positivity (%)'] <= group_75)][key]

        return list(middle_group_samples), group_25, group_75

    def plotVarVersusMean(self):
        dfFileName = sys.argv[1]
        plotName = sys.argv[2]
        crispGenesFileName = sys.argv[3]

        #####################
        with open(dfFileName, 'r') as f:
            df = pd.read_csv(f, header=0, sep=',')
        f.close()

        if 'env' in df.columns:
            df.drop(columns=['env'], inplace=True)

        if 'oro_thresh' in df.columns:
            df.drop(columns=['oro_thresh'], inplace=True)

        means_dict = dict(df.mean())
        vars_dict = dict(df.var())

        means_list = list(means_dict.values())
        vars_list = list(vars_dict.values())
        #####################

        #####################
        with open(crispGenesFileName, 'r') as f:
            crispGenes = f.read().splitlines()
        f.close()
        print(crispGenes)
        #####################

        #####################
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_xlabel('mean', fontsize=18)
        ax.set_ylabel('variance', fontsize=18)
        fig.suptitle(plotName, fontsize=24)

        plt.scatter(x=np.log10(means_list), y=np.log10(vars_list))

        crisp_x = list()
        crisp_y = list()

        for gene in crispGenes:
            crisp_x.append(np.log10(means_dict[gene]))
            crisp_y.append(np.log10(vars_dict[gene]))

        plt.scatter(x=crisp_x, y=crisp_y, marker='*', color='red')

        plt.savefig(plotName + '.png')
        #####################

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
                expr_row = expr_df[expr_df['sample'] == sample].drop(columns=['sample'])
                noise = np.random.normal(0, var, expr_row.shape)
                noised_expr_row = expr_row + noise
                noised_expr_row[noised_expr_row < 0 ] = 0
                myRand = random.randint(0, 1000000)
                while myRand in myRands:
                    myRand = random.randint(0, 1000000)
                myRands.append(myRand)
                new_sample = sample + '_' + str(threadID) + '_' + str(myRand)
                self.samples.append(new_sample)
                #new_sample = sample + '_' + str((i+1) * n * threadID)

                noised_expr_row['sample'] = new_sample
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

    def amplify_expr(self, df=None, metaDF = None, amplify = None, key='sample', numProcs=1):
        if amplify['n'] == 0:
            return df, metaDF

        expr_df = df
        meta_df = metaDF
        random.seed(amplify['seed'])

        # amplify set of samples that match column value
        q = Queue()
        processList = list()
        for i in range(numProcs):
            p = Process(target=self.process_samples_subset, args=(q,list(expr_df['sample']), expr_df['sample'], expr_df,
                                                                  meta_df, amplify['n'], amplify['var'], i, numProcs, amplify['seed'], key, ))
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

        #self.save_meta(meta_df, self.outputDir + '/metadata_amplify-' + str(amplify['n']) + '-' + str(amplify['var']) + '.csv')
        #self.save_expr(expr_df, self.outputDir + '/' + self.expr_outputfile_prefix + '_amplify-' + str(amplify['n']) + '-' + str(amplify['var']) + '.csv', transpose=True, cur_index_col='sample', new_index_col='gene')
        self.expr_outputfile_prefix += '_amplify-' + str(amplify['n']) + '-' + str(amplify['var']) + '_'

        return expr_df, meta_df

if __name__ == "__main__":
    main()
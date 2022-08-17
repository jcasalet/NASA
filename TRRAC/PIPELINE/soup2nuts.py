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


pd.options.mode.chained_assignment = None  # default='warn'

env_list= ['study', 'dissection', 'libprep']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expr', help='expression file list', default=None, required=True)
    parser.add_argument('-i', '--idir', help='input dir', default=None, required=True)
    parser.add_argument('-o', '--odir', help='output dir', default=None, required=True)
    parser.add_argument('-m', '--meta', help='metadata file name', default=None, required=True)
    args = parser.parse_args()

    # initialize expression and meta data from files
    myrnaseqdata = RNASeqData(rnaSeqFileList=args.expr, metaFileList=args.meta, inputDir=args.idir, outputDir=args.odir,
                              normalize = 'afterMerge', standardize = 'never', log_xform=None, oro_thresholds_per_study=True)

    # convert gene ids to gene names
    print('dims before converting to names: ', myrnaseqdata.expressionDF.shape)
    myrnaseqdata.convertIdsToNames()
    print('dims after converting to names: ', myrnaseqdata.expressionDF.shape)

    # filter only protein-coding genes
    print('dims before filter by type: ', myrnaseqdata.expressionDF.shape)
    myrnaseqdata.filterGenesByType(gene_type='protein_coding', id='gene')
    print('dims after filter by type: ', myrnaseqdata.expressionDF.shape)

    # filter genes with 0 count in at least n samples
    print('dims before filter 0: ', myrnaseqdata.expressionDF.shape)
    myrnaseqdata.filterGenesByPercentZeroCount(n=80)
    print('dims after filter 0: ', myrnaseqdata.expressionDF.shape)

    # reduce number of genes to n
    print('dims before filter by top n: ', myrnaseqdata.expressionDF.shape)
    myrnaseqdata.filterGenesByTopNSD(n=0)
    print('dims after filter by top n: ', myrnaseqdata.expressionDF.shape)

    # amplify number of samples by n more samples
    print('dims before amplify: ', myrnaseqdata.expressionDF.shape)
    myrnaseqdata.amplify_expr(n=0, var=50, seed=23, key='sample', numProcs=4)
    print('dims after amplify: ', myrnaseqdata.expressionDF.shape)

    # vst can only be run on count data, NOT on zscores (and need to drop any meta cols from expr like env, oro_thresh
    # use DESeq to shrink data with vsd, rld, and ntd transformations
    myrnaseqdata.shrinkExpression('vsd', myrnaseqdata.outputDir + '/expr_vsd.csv')

    # prepare data for crisp consumption
    myrnaseqdata.prep4Crisp(inputFile=None,
                            outputFile=myrnaseqdata.outputDir + '/' + myrnaseqdata.outputFilePrefix + '.pkl',
                            env_list=env_list, target='oro_thresh')

    myrnaseqdata.prep4Crisp(inputFile=myrnaseqdata.outputDir + '/expr_vsd.csv',
                            outputFile=myrnaseqdata.outputDir + '/' + myrnaseqdata.outputFilePrefix + '_vsd.pkl',
                            env_list = env_list, target='oro_thresh')

class RNASeqData():
    def __init__(self, inputDir ='.',
                 outputDir='.',
                 outputFilePrefix='expr',
                 rnaSeqFileList=None,
                 metaFileList=None,
                 oro_scale='raw',
                 oro_thresholds_per_study=True,
                 middle50_samples=True,
                 RScriptPath='/usr/local/bin/Rscript',
                 normalize = 'afterMerge',
                 standardize='never',
                 log_xform=None):
        self.RScriptPath = RScriptPath
        self.inputDir = inputDir
        self.outputDir = outputDir
        self.outputFilePrefix = outputFilePrefix
        #self.metaFileName = metaFileName
        self.normalize = normalize
        self.standardize = standardize
        self.log_xform = log_xform
        self.oro_thresholds_per_study=oro_thresholds_per_study
        self.oro_scale=oro_scale
        self.env_list = env_list

        # set up metadata
        self.metaFileList = metaFileList.split(',')
        self.metaDict = dict()
        for f in self.metaFileList:
            self.metaDict[f] = pd.read_csv(self.inputDir + '/' + f, sep=',', header=0)
        self.metaDF = pd.concat(list(self.metaDict.values()), ignore_index=True)
        self.metaDF.columns = list(map(lambda i: i.lower(), self.metaDF.columns))

        # set up expression df dict
        self.rnaSeqFileList = rnaSeqFileList.split(',')
        self.rnaExprDataDict = dict()
        if self.normalize == 'afterMerge':
            for f in self.rnaSeqFileList:
                self.rnaExprDataDict[f] = pd.read_csv(self.inputDir + '/' + f, sep=',', header=0)
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            inputFile=self.outputDir + '/expr_before_normalize.csv'
            outputFile=self.outputDir + '/expr_after_normalize.csv'
            metaFile=self.outputDir + '/meta.csv'
            self.save_expr(inputFile)
            self.save_meta(metaFile)
            cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/normalize.R', '--args', inputFile, metaFile, outputFile]
            self.callR(cmd)
            self.expressionDF = pd.read_csv(outputFile, sep=',', header=0)
            self.outputFilePrefix += '_norm-after-merge_'

        elif self.normalize == 'beforeMerge' and self.standardize == 'beforeMerge':
            for e, m in zip(self.rnaSeqFileList, self.metaFileList):
                cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/normalize.R', '--args', self.inputDir + '/' + e,
                       self.inputDir + '/' + m, self.outputDir + '/mor_' + e]
                self.callR(cmd)
                self.rnaExprDataDict[e] = pd.read_csv(self.outputDir + '/mor_' + e, sep=',', header=0)
                # there's an issue in R where it writes  either  an unnamed column for genes or 2 columns: one named 'gene' and one named 'Unnamed: 0'
                if 'gene' in self.rnaExprDataDict[e].columns and 'Unnamed: 0' in self.rnaExprDataDict[e].columns:
                    self.rnaExprDataDict[e].drop(columns=['Unnamed: 0'], inplace=True)
            for e in self.rnaExprDataDict.keys():
                genes=list(self.rnaExprDataDict[e]['gene'])
                self.rnaExprDataDict[e].drop(columns=['gene'], inplace=True)
                numeric_cols = list(self.rnaExprDataDict[e].columns)
                self.rnaExprDataDict[e] = self.rnaExprDataDict[e].apply(zscore)
                self.rnaExprDataDict[e]['gene'] = genes
                self.rnaExprDataDict[e] = self.rnaExprDataDict[e][['gene'] + numeric_cols]
                #self.rnaExprDataDict[e] = (self.rnaExprDataDict[e] - self.rnaExprDataDict[e].mean(axis=1)) / self.rnaExprDataDict[e].std(axis=1)

            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            self.outputFilePrefix += '_norm-before-merge_std-before-merge_'

        elif self.normalize == 'beforeMerge':
            for e, m in zip(self.rnaSeqFileList, self.metaFileList):
                cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/normalize.R', '--args', self.inputDir + '/' + e,
                       self.inputDir + '/' + m, self.outputDir + '/mor_' + e]
                self.callR(cmd)
                self.rnaExprDataDict[e] = pd.read_csv(self.outputDir + '/mor_' + e, sep=',', header=0)
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            self.outputFilePrefix += '_norm-before-merge_'

        elif self.normalize == 'never':
            for f in self.rnaSeqFileList:
                self.rnaExprDataDict[f] = pd.read_csv(self.inputDir + '/' + f, sep=',', header=0)
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), list(self.rnaExprDataDict.values()))
            self.outputFilePrefix += '_unnorm_'

        else:
            print('no such normalize value: ', self.normalize)
            sys.exit(1)

        self.genes = list(self.expressionDF['gene'])

        if not self.log_xform is None:
            self.my_log(base=self.log_xform)

        # set up samples
        self.samples = list(self.metaDF['sample'])
        self.oro_thresholds_per_study_dict=dict()
        if middle50_samples or self.oro_thresholds_per_study:
            for group in set(self.metaDF['group']):
                for study in set(self.metaDF['study']):
                    m50,u25,l75 = self.getMiddle50('sample', group, study)
                    self.oro_thresholds_per_study_dict[study] = (u25 + l75) / 2
                    self.samples += m50

        # subset expressionDF based on sample list
        self.expressionDF = self.expressionDF[self.expressionDF.columns.intersection(self.samples)]

        # add back the gene column
        self.expressionDF['gene'] = self.genes
        cols = list(self.expressionDF.columns)
        cols = [cols[-1]] + cols[:-1]
        self.expressionDF = self.expressionDF[cols]

        # permute samples to be in same order in expr as in meta
        self.permuteSamples()

        # convert genes x samples to samples x genes
        self.expressionDF = self.transpose_df(self.expressionDF, 'gene', 'sample')

        # combine same gene ids into one row
        print('dims before collapse: ', self.expressionDF.shape)
        self.collapseGeneCounts()
        print('dims after collapse: ', self.expressionDF.shape)


    def my_log(self, base=10):
        genes = list(self.expressionDF['gene'])
        df_np = np.array(self.expressionDF.drop(columns=['gene']))
        if base == 10:
            df_np_log = pd.DataFrame(np.log10(df_np+1))
        elif base == 2:
            df_np_log = pd.DataFrame(np.log2(df_np+1))
        else:
            print('base unknown: ', base)
            sys.exit(1)
        df_np_log['gene'] = genes
        cols = list(self.expressionDF.columns[1:])
        df_np_log.columns = cols + ['gene']
        df_np_log = df_np_log[['gene'] + cols]
        self.expressionDF = df_np_log
        self.outputFilePrefix += '_log' + str(base) + '_'


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

    def prep4Crisp(self, inputFile, outputFile, env_list, target):
        if not inputFile is None:
            self.expressionDF = self.read_expr(inputFile)
        self.add_env(env_list=env_list)
        if target == 'oro_thresh':
            self.add_oro()
        else:
            self.setTargetByKey(target)
        self.save_expr(outputFile)

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

    def save_expr(self, fileName, transpose=False, dropCols=[], cur_index_col=None, new_index_col=None):
        if transpose:
            df = self.transpose_df(self.expressionDF.drop(columns=dropCols), cur_index_col, new_index_col)
            if fileName.endswith('.csv'):
                df.to_csv(fileName, sep=',', index=None)
            elif fileName.endswith('.pkl'):
                df.to_pickle(fileName)
            else:
                print('unknown filename extension: ', fileName)
                sys.exit(1)
        else:
            df = self.expressionDF.drop(columns=dropCols)
            if fileName.endswith('.csv'):
                df.to_csv(fileName, sep=',', index=None)
            elif fileName.endswith('.pkl'):
                df.to_pickle(fileName)
            else:
                print('unknown filename extension: ', fileName)
                sys.exit(1)

    def save_meta(self, fileName):
        self.metaDF.to_csv(fileName, sep=',', index=None)

    def convertIdsToNames(self):
        input_to_R = self.outputDir + '/expr_input.csv'
        output_from_R = self.outputDir + '/expr_output.csv'
        self.save_expr(input_to_R, transpose=True, dropCols=[], cur_index_col='sample', new_index_col='gene')
        R_cmd = ['/usr/local/bin/R', '-f', '/Users/jcasalet/convert_id_to_gene.R', '--args', input_to_R, output_from_R]
        self.callR(R_cmd)
        self.expressionDF = self.read_expr(output_from_R)

    def shrinkExpression(self, transformation, fileName):
        if not transformation in ['vsd', 'rld', 'ntd']:
            print('transformation: ' + transformation + ' not known')
            sys.exit(1)
        input_to_R = self.outputDir + '/expr_before_shrink.csv'
        output_from_R = fileName
        meta_file = self.outputDir + '/meta.csv'
        self.save_expr(input_to_R, transpose=True, dropCols=[], cur_index_col='sample', new_index_col='gene')
        #self.save_expr(input_to_R)
        self.save_meta(meta_file)
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

    def my_standardize(self, by):
        self.outputFilePrefix += '_std_'
        #self.expressionDF = (self.expressionDF - self.expressionDF.mean())/self.expressionDF.std()
        if by == 'col':
            self.expressionDF = (self.expressionDF - self.expressionDF.mean(axis=1)) / self.expressionDF.std(axis=1)
        elif by == 'row':
            self.expressionDF = (self.expressionDF - self.expressionDF.mean(axis=0)) / self.expressionDF.std(axis=0)

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

    def add_env(self, env_list=None):
        # create dictionary with sample id as key and concatenated column strings as value
        env_dict = dict()
        for i in range(len(self.metaDF)):
            key = self.metaDF.iloc[i]['sample']
            if not env_list:
                value = self.metaDF.iloc[i]['study'] + ':' + self.metaDF.iloc[i]['dissection'] + ':' + self.metaDF.iloc[i]['libprep']
            else:
                counter = 0
                for e in env_list:
                    if counter == 0:
                        value = self.metaDF.iloc[i][e]
                        counter+=1
                    else:
                        value = value + ':' + self.metaDF.iloc[i][e]
            env_dict[key] = value

        # join env dictionary to data frame
        self.expressionDF['env'] = self.expressionDF['sample'].map(env_dict)

    def filterGenesByPercentZeroCount(self, n=0):
        if n == 0:
            pass
        else:
            df = self.transpose_df(self.expressionDF, 'sample', 'gene')
            df = df[(df == 0).sum(axis='columns') <= int(n * len(df.columns))]
            self.expressionDF = self.transpose_df(df, 'gene', 'sample')
            self.genes = list(self.expressionDF.columns)[1:]

    def filterGenesByType(self, gene_type='protein_coding', id='id'):
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
            self.expressionDF = self.expressionDF[self.expressionDF.columns.intersection(filter_columns)]
            new_columns = list(self.expressionDF.drop(columns=['sample']))
            #gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene name'])
            gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene stable ID'])
            self.expressionDF.columns = ['sample'] + gene_names
        elif id == 'gene':
            filter_genes = list(gene_info[gene_info['Gene type'] == gene_type]['Gene name'])
            filter_columns = ['sample'] + filter_genes
            self.expressionDF = self.expressionDF[self.expressionDF.columns.intersection(filter_columns)]
        self.expressionDF = self.expressionDF.loc[:, self.expressionDF.columns.notna()]
        self.genes = list(self.expressionDF.columns)[1:]
        self.outputFilePrefix += '_protein-coding_'

    def filterGenesByTopNSD(self, n):
        # df is genes X samples
        # calculate var, sort cols into n highest vars, drop shape[1]-n cols
        # first find range of var and print to stdout
        if n == 0:
            pass
        else:
            df = self.transpose_df(self.expressionDF, 'sample', 'gene')
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
            self.expressionDF = self.transpose_df(df, 'gene', 'sample')
        self.genes = list(self.expressionDF.columns)[1:]
        self.outputFilePrefix += '_top-' + str(n) + '_'

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
        self.genes = list(self.expressionDF.columns)[1:]


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

    def amplify_expr(self, n, var, seed=0, key='sample', numProcs=1):
        if n == 0:
            return

        expr_df = self.expressionDF
        meta_df = self.metaDF
        random.seed(seed)

        # amplify set of samples that match column value
        q = Queue()
        processList = list()
        for i in range(numProcs):
            p = Process(target=self.process_samples_subset, args=(q,list(expr_df['sample']), expr_df['sample'], expr_df, meta_df, n, var, i, numProcs, seed, key, ))
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

        '''genes = expr_df['gene']
        expr_df = expr_df.drop(columns=['gene'])
        expr_df = np.clip(expr_df, 0, a_max=None)
        expr_df.insert(0, 'gene', genes)'''
        self.expressionDF = expr_df
        self.metaDF = meta_df
        self.outputFilePrefix += '_amplify-' + str(n) + '-' + str(var) + '_'

if __name__ == "__main__":
    main()
import subprocess
import pandas as pd
from statistics import mean, median
from sklearn import preprocessing
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
from pybiomart import Server

pd.options.mode.chained_assignment = None  # default='warn'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expr', help='expression file list', default=None, required=True)
    parser.add_argument('-i', '--idir', help='input dir', default=None, required=True)
    parser.add_argument('-o', '--odir', help='output dir', default=None, required=True)
    parser.add_argument('-m', '--meta', help='metadata file name', default=None, required=True)
    args = parser.parse_args()
    myrnaseqdata = RNASeqData(rnaSeqFileList=args.expr, metaFileName=args.meta, inputDir=args.idir, outputDir=args.odir,
                              normalizeFirst = False, oro_thresholds_per_study=True, env_list=['study', 'dissection', 'libprep'])
    myrnaseqdata.filterGenesByType(gene_type='protein_coding')
    myrnaseqdata.collapseGeneCounts()
    myrnaseqdata.filterGenesByCount(0.7)
    myrnaseqdata.add_env(env_list=['group'])
    myrnaseqdata.add_oro()
    myrnaseqdata.save_expr('expr.csv', transpose=True)
    myrnaseqdata.save_meta('meta.csv', transpose=False)
    myrnaseqdata.callR('/usr/local/bin/R', '/Users/jcasalet/Desktop/my_test.R', myrnaseqdata.outputDir + '/expr.csv',
                       myrnaseqdata.outputDir + '/meta.csv')

    print(myrnaseqdata.samples)
    print(myrnaseqdata.oro_thresholds_per_study_dict)
    print(myrnaseqdata.expressionDF['oro_thresh'])
    print(myrnaseqdata.expressionDF['env'])

class RNASeqData():
    def __init__(self, inputDir ='.', rnaSeqFileList=None, metaFileName=None, outputDir='.',
                 oro_scale='raw', middle50_samples=False, RScriptPath='/usr/local/bin/Rscript', normalizeFirst = False,
                 oro_thresholds_per_study=True, env_list=None):
        self.RScriptPath = RScriptPath
        self.inputDir = inputDir
        self.outputDir = outputDir
        self.metaFileName = metaFileName
        self.oro_scale=oro_scale
        self.normalizeFirst = normalizeFirst
        self.oro_thresholds_per_study=oro_thresholds_per_study
        self.env_list = env_list

        # set up expression df
        self.rnaSeqFileList = rnaSeqFileList.split(',')
        dfList = list()
        for f in self.rnaSeqFileList:
            dfList.append(pd.read_csv(self.inputDir + '/' + f, sep=',', header=0))
        import functools as ft
        if not normalizeFirst:
            self.expressionDF = ft.reduce(lambda left, right: pd.merge(left, right, on='gene'), dfList)

        self.genes = list(self.expressionDF['gene'])

        # set up metadata
        self.metaDF = pd.read_csv(inputDir + '/' + self.metaFileName, sep=',', header=0)
        self.metaDF.columns = list(map(lambda i: i.lower(), self.metaDF.columns))

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

        # convert genes x samples to samples x genes
        self.expressionDF = self.transpose_df(self.expressionDF, 'gene', 'sample')

    def save_expr(self, fileName, transpose):
        if transpose:
            df = self.transpose_df(self.expressionDF.drop(columns=['env', 'oro_thresh']), 'sample', 'gene')
            df.to_csv(self.outputDir +  '/' + fileName, sep=',', index=None)
        else:
            self.expressionDF.drop(columns=['env', 'oro_thresh']).to_csv(self.outputDir +  '/' + fileName, sep=',', index=None)

    def save_meta(self, fileName, transpose):
        self.metaDF.to_csv(self.outputDir + '/' + fileName, sep=',', index=None)

    def callR(self, RBinaryPath, RScriptPath, exprFile, metaFile):
        import subprocess
        cmd = [RBinaryPath, '-f', RScriptPath, '--args', exprFile, metaFile, '--no-save']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        o, e = proc.communicate(timeout=60)

        print('Output: ' + o.decode('ascii'))
        #print('Error: ' + e.decode('ascii'))
        print('code: ' + str(proc.returncode))

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
        '''rownames = list(df[cur_index_col])
        colnames = list(df.drop(columns=[cur_index_col]).columns)
        df.drop(columns = [cur_index_col], inplace=True)
        df=df.transpose()
        df.columns = rownames
        df[new_index_col] = colnames
        cols = [new_index_col] + list(df.columns[:-1])
        df = df[cols]
        df.reset_index(inplace=True)'''
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

    def filterGenesByCount(self, alpha):
        print('dims before filter 0: ', self.expressionDF.shape)
        df = self.transpose_df(self.expressionDF, 'sample', 'gene')
        df = df[(df == 0).sum(axis='columns') <= int(alpha * len(df.columns))]
        self.expressionDF = self.transpose_df(df, 'gene', 'sample')
        print('dims after filter 0: ', self.expressionDF.shape)

    def filterGenesByType(self, gene_type='protein_coding'):
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
        filter_genes = list(gene_info[gene_info['Gene type'] == gene_type]['Gene stable ID'])
        filter_columns = ['sample'] + filter_genes
        self.expressionDF = self.expressionDF[self.expressionDF.columns.intersection(filter_columns)]
        new_columns = list(self.expressionDF.drop(columns=['sample']))
        gene_names = list(gene_info[gene_info['Gene stable ID'].isin(new_columns)]['Gene name'])
        self.expressionDF.columns = ['sample'] + gene_names
        self.expressionDF = self.expressionDF.loc[:, self.expressionDF.columns.notna()]

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
            df = df.append(collapsed_genes[gene], ignore_index=True)
            index=len(self.expressionDF)-1
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

if __name__ == "__main__":
    main()
import numpy as np
from sklearn.svm import SVR
#from sklearn.model_selection import train_test_split
import pandas as pd
from sklearn.model_selection import cross_val_score
from statistics import mean
import math
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import scale
from pandas import DataFrame
import pickle
from sklearn.ensemble import RandomForestRegressor
#from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
import time
from scipy import stats
import argparse
from sklearn.linear_model import ElasticNet
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--chr", action="store", dest="chr",
                    help="specify the chromosome number")
parser.add_argument("--training_pop", action="store", dest="training_pop",
                    help="Imputation training population name")
parser.add_argument("--testing_pop", action="store", dest="testing_pop",
                    help="Imputation testing population name")
parser.add_argument("--output_dir", action="store", dest="output_dir",
                    help="specify the output directory. Start and end with slash")
parser.add_argument("--trn_data_path", action="store", dest="trn_data_path",
                    help="Specify train data path. Start and end with slash")
parser.add_argument("--tst_data_path", action="store", dest="tst_data_path",
                    help="Specify test data path. Start and end with slash")
parser.add_argument("--param_file_path", action="store", dest="param_file_path",
                    help="Specify path to the gridsearch hyperparameter files")

args = parser.parse_args()
chrom = str(args.chr)
trn_pop = str(args.training_pop)
tst_pop = str(args.testing_pop)
output = str(args.output_dir)
trn_data_path = str(args.trn_data_path)
tst_data_path = str(args.tst_data_path)
param_file_path = str(args.param_file_path)

#important functions needed
def get_filtered_snp_annot (snpfilepath):
     snpanot = pd.read_csv(snpfilepath, sep="\t")
     snpanot = snpanot[(((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) |
                        ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) |
                        ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) |
                        ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) |
                        ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) |
                        ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) |
                        ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) |
                        ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))) &
                       (snpanot["rsid"].notna())]
     snpanot = snpanot.drop_duplicates(["varID"])
     return snpanot


def get_gene_annotation (gene_anot_filepath, chrom, gene_types=["protein_coding"]):
     gene_anot = pd.read_csv(gene_anot_filepath, sep="\t")
     gene_anot = gene_anot[(gene_anot["chr"]==str(chrom)) &
                           (gene_anot["gene_type"].isin(gene_types))]
     return gene_anot

def get_gene_type (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_type = gene_type.iloc[0,5]
     return gene_type

def get_gene_name (gene_anot, gene):
     gene_name = gene_anot[gene_anot["gene_id"]==gene]
     gene_name = gene_name.iloc[0,2]
     return gene_name

def get_gene_coords (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_coord = [gene_type.iloc[0,3], gene_type.iloc[0,4]]
     return gene_coord

def get_covariates (cov_filepath):
      cov = pd.read_csv(cov_filepath, sep=" ")
      cov = cov.set_index("IID") #make IID to be the row names
      cov.index.names = [None] # remove the iid name from the row
      pc = ["PC1", "PC2", "PC3"] #a list of the PCs to retain
      cov = cov[pc]
      return cov

def get_gene_expression(gene_expression_file_name, gene_annot):
	expr_df = pd.read_csv(gene_expression_file_name, header = 0, index_col = 0, delimiter='\t')
	expr_df = expr_df.T
	inter = list(set(gene_annot['gene_id']).intersection(set(expr_df.columns)))
	#print(len(inter))
	expr_df = expr_df.loc[:, inter ]
	return expr_df

def adjust_for_covariates (expr_vec, cov_df):   
      reg = LinearRegression().fit(cov_df, expr_vec)
      ypred = reg.predict(cov_df)
      residuals = expr_vec - ypred
      residuals = scale(residuals)
      return residuals

def get_maf_filtered_genotype(genotype_file_name,  maf=0.01):
	gt_df = pd.read_csv(genotype_file_name, 'r', header = 0, index_col = 0,delimiter='\t')
	effect_allele_freqs = gt_df.mean(axis=1)
	effect_allele_freqs = [ x / 2 for x in effect_allele_freqs ]
	effect_allele_boolean = pd.Series([ ((x >= maf) & (x <= (1 - maf))) for x in effect_allele_freqs ]).values
	gt_df = gt_df.loc[ effect_allele_boolean ]
	gt_df = gt_df.T
	return gt_df

def get_cis_genotype (gt_df, snp_annot, coords, cis_window=1000000):
      snp_info = snpannot[(snpannot['pos'] >= (coords[0] - cis_window)) &
                          (snpannot['rsid'].notna()) & (snpannot['pos'] <= (coords[1] + cis_window))]
      if len(snp_info) == 0:
          return 0
      else:
           gtdf_col = list(gt_df.columns)
           snpinfo_col = list(snp_info["varID"])
           intersect = snps_intersect(gtdf_col, snpinfo_col) #this function was defined earlier
           cis_gt = gt_df[intersect]
           return cis_gt


def snps_intersect(list1, list2):
     return list(set(list1) & set(list2))


#train data files
#AFA_chr1.gz
#AFA_GEX.gz
#AFA_chr1_annot.gz
#gencodev18.gz
#AFA_3_PC.txt
snp_dosage_file = trn_data_path+trn_pop.upper()+"_chr"+chrom+".txt"
gene_expression_file = trn_data_path+trn_pop.upper()+"_GEX.txt"
pc_file = trn_data_path+trn_pop.upper()+"_3_PCs.txt"
gene_annotation_file = trn_data_path+"/gencode.v18.annotation.parsed.txt"
snp_annotation_file = trn_data_path+trn_pop.upper()+"_chr"+chrom+"_annot.txt"


#test data files
test_snp_dosage_file = tst_data_path+tst_pop.upper()+"_chr"+chrom+".txt"
test_gene_expression_file = tst_data_path+tst_pop.upper()+"_GEX.txt"
test_pc_file = tst_data_path+tst_pop.upper()+"_3_PCs.txt"
test_gene_annotation_file = tst_data_path+"/gencode.v28.annotation.parsed.txt"
test_snp_annotation_file = tst_data_path+tst_pop.upper()+"_chr"+chrom+"_annot.txt"

#train functioning
snpannot = get_filtered_snp_annot(snp_annotation_file)
geneannot = get_gene_annotation(gene_annotation_file, chrom)
expr_df = get_gene_expression(gene_expression_file, geneannot) #this had to created early to avoid empty df downstream
annot_geneid = geneannot["gene_id"]#remove decimal from gene_id
annot_geneid = list(annot_geneid)
agid = []
for i in annot_geneid:
	agid.append(i[0:(i.find("."))])
geneannot["gene_id"] = agid #replace with non decimal gene_id

cov = get_covariates(pc_file)

genes = list(expr_df.columns)
gt_df = get_maf_filtered_genotype(snp_dosage_file)
train_ids = list(gt_df.index)
train_g = [] #where to store the non decimal gene_id
for i in genes:
	train_g.append(i[0:(i.find("."))])
expr_df.columns = train_g #use the non decimal gene_id to rename the expr_df col
genes = list(expr_df.columns) #take out the new non decimal gene_id
#adj_exp_frame = pd.DataFrame()


#test functioning
test_snpannot = get_filtered_snp_annot(test_snp_annotation_file)
test_geneannot = get_gene_annotation(test_gene_annotation_file, chrom)
test_expr_df = get_gene_expression(test_gene_expression_file, test_geneannot) #important to create early
test_annot_geneid = test_geneannot["gene_id"]#remove decimal from gene_id
test_annot_geneid = list(test_annot_geneid)
test_agid = []
for i in test_annot_geneid:
	test_agid.append(i[0:(i.find("."))])
test_geneannot["gene_id"] = test_agid #replace with non decimal gene_id

test_cov = get_covariates(test_pc_file)

test_genes = list(test_expr_df.columns)
test_gt_df = get_maf_filtered_genotype(test_snp_dosage_file)
test_ids = list(test_gt_df.index)
test_g = []
for i in test_genes:
	test_g.append(i[0:(i.find("."))])
test_expr_df.columns = test_g
test_genes = list(test_expr_df.columns)

#frame to store the ypred and test adjusted expression
ypred_frame_rf = pd.DataFrame()
ypred_frame_svr = pd.DataFrame()
ypred_frame_knn = pd.DataFrame()


#read in the grid search best result files and take the params to fit the model
rf_grid = pd.read_csv(param_file_path+trn_pop.upper()+"_best_grid_rf_chr"+chrom+
                      "_full.txt", sep="\t")

knn_grid = pd.read_csv(param_file_path+trn_pop.upper()+"_best_grid_knn_chr"+chrom+
                       "_full.txt", sep="\t")

svr_grid = pd.read_csv(param_file_path+trn_pop.upper()+"_best_grid_svr_chr"+chrom+
                       "_full.txt", sep="\t")


#text file where to write out the cv and test results
open(output+"grid_optimized_"+trn_pop.upper()+"_2_"+tst_pop.upper()+"_rf_cor_test_chr"+str(chrom)+
     ".txt", "w").write("gene_id"+"\t"+"gene_name"+"\t"+"pearson_yadj_vs_ypred (a)"+"\t"+"a_pval"+
                        "\t"+"pearson_yobs_vs_ypred (b)"+"\t"+"b_pval"+"\t"+"spearman_yadj_vs_ypred (c)"+
                        "\t"+"c_pval"+"\t"+"spearman_yobs_vs_ypred (d)"+"\t"+"d_pval"+"\n")

open(output+"grid_optimized_"+trn_pop.upper()+"_2_"+tst_pop.upper()+"_knn_cor_test_chr"+str(chrom)+
     ".txt", "w").write("gene_id"+"\t"+"gene_name"+"\t"+"pearson_yadj_vs_ypred (a)"+"\t"+"a_pval"+"\t"+
                        "pearson_yobs_vs_ypred (b)"+"\t"+"b_pval"+"\t"+"spearman_yadj_vs_ypred (c)"+"\t"+"c_pval"+"\t"+
                        "spearman_yobs_vs_ypred (d)"+"\t"+"d_pval"+"\n")

open(output+"grid_optimized_"+trn_pop.upper()+"_2_"+tst_pop.upper()+"_svr_cor_test_chr"+str(chrom)+
     ".txt", "w").write("gene_id"+"\t"+"gene_name"+"\t"+"pearson_yadj_vs_ypred (a)"+"\t"+"a_pval"+"\t"+
                        "pearson_yobs_vs_ypred (b)"+"\t"+"b_pval"+"\t"+"spearman_yadj_vs_ypred (c)"+"\t"+"c_pval"+"\t"+
                        "spearman_yobs_vs_ypred (d)"+"\t"+"d_pval"+"\n")

for gene in genes:
     
    if gene in test_genes:
        
        train_coords = get_gene_coords(geneannot, gene)        
        test_coords = get_gene_coords(test_geneannot, gene)
        gene_name = get_gene_name(geneannot, gene)
        
        
        expr_vec = expr_df[gene]#observed exp
        test_expr_vec = test_expr_df[gene]#observed exp
        
        adj_exp = adjust_for_covariates(list(expr_vec), cov)#adjusted exp
        test_adj_exp = adjust_for_covariates(list(test_expr_vec), test_cov)#adjusted exp
        
        
        cis_gt = get_cis_genotype(gt_df, snpannot, train_coords)
        test_cis_gt = get_cis_genotype(test_gt_df, test_snpannot, test_coords)

        if (type(cis_gt) != int) & (type(test_cis_gt) != int):#just to be sure the cis genotype is not empty
             
             gg = [gene] #just to cast the gene id to list because pandas need it to be in list before it can be used as col name

             #take the snps
             train_snps = list(cis_gt.columns)
             test_snps = list(test_cis_gt.columns)
             snp_intersect = snps_intersect(train_snps, test_snps)

             cis_gt = cis_gt[snp_intersect]
             test_cis_gt = test_cis_gt[snp_intersect]

             if (cis_gt.shape[1] > 0) & (test_cis_gt.shape[1] > 0): #make sure that the cis_gt is not empty
                  
             
                  #build the model
                  #adj_exp = adj_exp.values #not needed after making adj_exp a numpy array above
                  cis_gt = cis_gt.values
                  test_cis_gt = test_cis_gt.values
                  test_yobs = test_expr_vec.values
                  
                  #Random Forest
                  gnlist = list(rf_grid['Gene_Name'])
                  f = gene_name in gnlist
                  if f != False: #Just to be sure that the gene exist in the RF best grid dataframe
                       
                       n_tree = rf_grid[rf_grid['Gene_Name']==gene_name].iloc[0,3]
                       rf = RandomForestRegressor(random_state=1234, n_estimators=n_tree)
                       rf.fit(cis_gt, adj_exp.ravel())
                       ypred = rf.predict(test_cis_gt)

                       #prepare ypred for writing out to a file
                       ypred_pd = pd.DataFrame(ypred)
                       
                       ypred_pd.columns = gg
                       ypred_pd.index = test_ids
                       ypred_frame_rf = pd.concat([ypred_frame_rf, ypred_pd], axis=1, sort=True)
                       
                       pa = stats.pearsonr(test_adj_exp, ypred)
                       pacoef = str(float(pa[0]))
                       papval = str(float(pa[1]))
                       pb = stats.pearsonr(test_yobs, ypred)
                       pbcoef = str(float(pb[0]))
                       pbpval = str(float(pb[1]))
                       sc = stats.spearmanr(test_adj_exp, ypred)
                       sccoef = str(float(sc[0]))
                       scpval = str(float(sc[1]))
                       sd = stats.spearmanr(test_yobs, ypred)
                       sdcoef = str(float(sd[0]))
                       sdpval = str(float(sd[1]))
                       open(output+"grid_optimized_"+trn_pop.upper()+"_2_"+tst_pop.upper()+
                            "_rf_cor_test_chr"+str(chrom)+".txt", "a").write(gene+"\t"+gene_name+"\t"+pacoef+"\t"
                                                                             +papval+"\t"+pbcoef+"\t"+pbpval+"\t"+sccoef+
                                                                             "\t"+scpval+"\t"+sdcoef+"\t"+sdpval+"\n")

                  #Support Vector Machine
                  gnlist = list(svr_grid['Gene_Name'])
                  f = gene_name in gnlist
                  if f != False: #Just to be sure that the gene exist in the SVR best grid dataframe
                       
                       kernel = svr_grid[svr_grid['Gene_Name']==gene_name].iloc[0,3]
                       degree = svr_grid[svr_grid['Gene_Name']==gene_name].iloc[0,4]
                       c = svr_grid[svr_grid['Gene_Name']==gene_name].iloc[0,5]
                       svr = SVR(gamma="scale", kernel=kernel, degree=degree, C=c)
                       svr.fit(cis_gt, adj_exp.ravel())
                       ypred = svr.predict(test_cis_gt)
                       
                       #prepare ypred for writing out to a file
                       yprep_pd = pd.DataFrame(ypred)
                       
                       ypred_pd.columns = gg
                       ypred_pd.index = test_ids
                       ypred_frame_svr = pd.concat([ypred_frame_svr, ypred_pd], axis=1, sort=True)
                       
                       pa = stats.pearsonr(test_adj_exp, ypred)
                       pacoef = str(float(pa[0]))
                       papval = str(float(pa[1]))
                       pb = stats.pearsonr(test_yobs, ypred)
                       pbcoef = str(float(pb[0]))
                       pbpval = str(float(pb[1]))
                       sc = stats.spearmanr(test_adj_exp, ypred)
                       sccoef = str(float(sc[0]))
                       scpval = str(float(sc[1]))
                       sd = stats.spearmanr(test_yobs, ypred)
                       sdcoef = str(float(sd[0]))
                       sdpval = str(float(sd[1]))
                       open(output+"grid_optimized_"+trn_pop.upper()+"_2_"+tst_pop.upper()+
                            "_svr_cor_test_chr"+str(chrom)+".txt", "a").write(gene+"\t"+gene_name+"\t"+pacoef+
                                                                              "\t"+papval+"\t"+pbcoef+"\t"+pbpval+"\t"+sccoef+
                                                                              "\t"+scpval+"\t"+sdcoef+"\t"+sdpval+"\n")
                  
                  
                  #K Nearest Neighbour
                  gnlist = list(knn_grid['Gene_Name'])
                  f = gene_name in gnlist
                  if f != False: #Just to be sure that the gene exist in the KNN best grid dataframe
                       
                       k = knn_grid[knn_grid['Gene_Name']==gene_name].iloc[0,3]
                       weight = knn_grid[knn_grid['Gene_Name']==gene_name].iloc[0,4]
                       knn = KNeighborsRegressor(n_neighbors=k, weights = weight)
                       knn.fit(cis_gt, adj_exp.ravel())
                       ypred = knn.predict(test_cis_gt)
                       
                       #prepare ypred for writing out to a file
                       yprep_pd = pd.DataFrame(ypred)
                       
                       ypred_pd.columns = gg
                       ypred_pd.index = test_ids
                       ypred_frame_knn = pd.concat([ypred_frame_knn, ypred_pd], axis=1, sort=True)
                       
                       pa = stats.pearsonr(test_adj_exp, ypred)
                       pacoef = str(float(pa[0]))
                       papval = str(float(pa[1]))
                       pb = stats.pearsonr(test_yobs, ypred)
                       pbcoef = str(float(pb[0]))
                       pbpval = str(float(pb[1]))
                       sc = stats.spearmanr(test_adj_exp, ypred)
                       sccoef = str(float(sc[0]))
                       scpval = str(float(sc[1]))
                       sd = stats.spearmanr(test_yobs, ypred)
                       sdcoef = str(float(sd[0]))
                       sdpval = str(float(sd[1]))
                       open(output+"grid_optimized_"+trn_pop.upper()+"_2_"+tst_pop.upper()+
                            "_knn_cor_test_chr"+str(chrom)+".txt", "a").write(gene+"\t"+gene_name+"\t"+pacoef+"\t"+papval+"\t"+pbcoef+
                                                                              "\t"+pbpval+"\t"+sccoef+"\t"+scpval+"\t"+sdcoef+
                                                                              "\t"+sdpval+"\n")

        

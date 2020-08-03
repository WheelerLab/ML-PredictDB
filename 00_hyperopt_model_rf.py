import numpy as np
from sklearn.svm import SVR
import pandas as pd
from sklearn.model_selection import cross_val_score
from statistics import mean
import math
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import scale
from pandas import DataFrame
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
import time
from scipy import stats
import argparse
from sklearn.linear_model import ElasticNet
#from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.metrics import explained_variance_score
from sklearn.metrics import r2_score
from sklearn.metrics import make_scorer#use to convert metrics to scoring callables
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
import hyperopt

r2 = make_scorer(r2_score, greater_is_better=True)


parser = argparse.ArgumentParser()
parser.add_argument("--chr", action="store", dest="chr",
                    help="specify the chromosome number")
parser.add_argument("--training_pop", action="store", dest="training_pop",
                    help="Imputation training population name")
parser.add_argument("--output_dir", action="store", dest="output_dir",
                    help="specify the output directory. Start and end with slash")
parser.add_argument("--data_path", action="store", dest="data_path",
                    help="Specify data path. Start and end with slash")
parser.add_argument("--evals", action="store", dest="chr",
                    help="specify the number of evals")

args = parser.parse_args()
chrom = str(args.chr)
pop = str(args.training_pop)
output = str(args.output_dir)
data_path = str(args.data_path)
max_evals= int(args.evals)

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
                          (snpannot['rsid'].notna()) &
                          (snpannot['pos'] <= (coords[1] + cis_window))]
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


#Set up hyperopt

#1 Define the objective function

def objective(params):
    regressor_type = params["type"]
    del params["type"]
    if regressor_type == "en":
        regressor = ElasticNet(max_iter=10000, random_state=1234, **params)
    elif regressor_type == "rf":
        regressor = RandomForestRegressor(random_state=1234, **params)
    elif regressor_type == "svm":
        regressor = SVR(gamma="scale", **params)
    elif regressor_type == "knn":
        regressor = KNeighborsRegressor(**params)
    else:
        return 0
    r2_mean = cross_val_score(regressor, x, y, scoring=r2, cv=5).mean()

    return {"loss": -r2_mean, "status": STATUS_OK}


#2 Define search space for each machine learning model

en_space = {
    "type": "en",
    "alpha": hp.lognormal("alpha", 1.0, 10.0)
}

rf_space = {
    "type": "rf",
    "n_estimators": hp.choice("trees", range(50, 550, 50))
}

svm_space = {
    "type": "svm",
    "C": hp.lognormal("C", 0, 1.0),
    "kernel": hp.choice("kernel", ["linear", "rbf", "sigmoid", "poly"]),
    "degree": hp.choice("degree", range(2,8,1))
}

knn_space = {
    "type": "knn",
    "n_neighbors": hp.choice("n_neighbors", range(3, 33, 2)),
    "weights": hp.choice("weights", ["uniform", "distance"]),
    "p": hp.choice("p", range(1, 4, 1))
}

#3 choose hyperopt search algorithm
algo = tpe.suggest #tpe = Tree-of-Parzen-Estimator
#chose max eval
#max_evals= 30

# Set file paths

#AFA_chr1.gz
#AFA_GEX.gz
#AFA_chr1_annot.gz
#gencodev18.gz
#AFA_3_PC.txt

snp_dosage_file = data_path+pop.upper()+"_chr"+chrom+".txt"
gene_expression_file = data_path+pop.upper()+"_GEX.txt"
pc_file = data_path+pop.upper()+"_3_PCs.txt"
gene_annotation_file = data_path+"/gencode.v18.annotation.parsed.txt"
snp_annotation_file = data_path+pop.upper()+"_chr"+chrom+"_annot.txt"

# parse the files

snpannot = get_filtered_snp_annot(snp_annotation_file)
geneannot = get_gene_annotation(gene_annotation_file, chrom)
cov = get_covariates(pc_file)
expr_df = get_gene_expression(gene_expression_file, geneannot)
genes = list(expr_df.columns)
gt_df = get_maf_filtered_genotype(snp_dosage_file)


# prepare files where to write the results for each ML method

#rf
open(output+pop+"_rf_hyperopt_chr"+chrom+
     ".txt", "w").write("gene_id"+"\t"+"gene_name"+"\t"+"chr"+"\t"+
                              "best_hyperparam"+"\t")

for i in range(0, max_evals, 1):
     open(output+pop+"_rf_hyperopt_chr"+chrom+
          ".txt", "a").write(str(i)+"\t")



#Go through all protein coding genes

for gene in genes:
    coords = get_gene_coords(geneannot, gene)
    gene_name = get_gene_name(geneannot, gene)
    expr_vec = expr_df[gene]
    
    adj_exp = adjust_for_covariates(list(expr_vec), cov)
    cis_gt = get_cis_genotype(gt_df, snpannot, coords)

    #build the model
    
    if (type(cis_gt) != int) & (cis_gt.shape[1] > 0):
         

         x = cis_gt.values
         y = adj_exp.ravel()


         #Random Forest
         trials = Trials() #reset the trials object
         best = fmin(fn=objective, space=rf_space, algo=algo,
                     max_evals=max_evals, trials=trials)
         result_table = pd.DataFrame(trials.results)
         best_hyperparam = hyperopt.space_eval(rf_space, best)
         best_hyperparam.pop("type") #just to remove "type" from the param dict
         
         open(output+pop+"_rf_hyperopt_chr"+
              chrom+".txt", "a").write("\n"+gene+"\t"+gene_name+"\t"+chrom+"\t"+str(best_hyperparam)+"\t")

         for i in range(0, max_evals, 1): #I negate the loss in order to get cvR2
              open(output+pop+"_rf_hyperopt_chr"+
                   chrom+".txt", "a").write(str(-1*result_table.loss[i])+"\t")

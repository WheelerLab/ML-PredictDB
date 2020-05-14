# ML-PredictDB
In this study we sought to optimize protein coding gene expression imputation performance across global populations in comparison to the popular transcriptome tool [PrediXcan](https://github.com/WheelerLab/PrediXcan). Specifically, we used non-linear machine learning (ML) models viz; random forest (RF), support vector regression (SVR), and K nearest neighbor (KNN) to build transcriptome imputation models, and evaluated their performance in comparison to elastic net (EN) - the algorithm used in [PrediXcan](https://github.com/WheelerLab/PrediXcan). We trained gene expression prediction models using genotype and blood monocyte transcriptome data from the Multi-Ethnic Study of Atherosclerosis (MESA) comprising individuals of African (AFA), Hispanic (HIS), and European (CAU) ancestries and tested them using genotype and whole blood transcriptome data from the Modeling the Epidemiology Transition Study (METS) comprising individuals of African ancestries. We show that the prediction performance is highest when the training and the testing population share similar ancestries regardless of the prediction algorithm used. While EN generally outperformed RF, SVR, and KNN, we found that RF outperforms EN for some genes, particularly between disparate ancestries, suggesting potential robustness and reduced variability of RF imputation performance across global populations. When applied to a high-density lipoprotein (HDL) phenotype, we show including RF prediction models in PrediXcan reveals potential gene associations missed by EN models. Therefore, by integrating non-linear modeling into PrediXcan and diversifying our training populations to include more global ancestries, we may uncover new genes associated with complex traits.

# Model Training and Testing
The model training and testing were done in two separate parts. Firstly, the elastic net (EN) models were trained and tested as described in [Mogil et al., 2018](https://github.com/WheelerLab/DivPop). Secondly and lastly, since the non-linear machine learning (ML) models do not have weights/beta values as we would have in a traditional regularization linear regression models like EN, which would have allowed us to build a DB file containing the weight values of the eQTLs in each gene model and do model testing whenever, we were limited to training and testing at the same time and in the same script. Following that analogy, the ML models were built as described below:
  1. We conducted the first training/testing using the training population via five-fold cross-validation. Specifically, we used the 00_gridsearch_model.py script to carry out a gridsearch through each ML algorithm hyperparameters in order to determine the hyperparameter values or combination of values that yields the optimum five-fold cross-validated imputation performance (R<sup>2</sup>) for each gene. The files containing the optimum hyperparameter values and (R<sup>2</sup>) for each gene across training population are KNN_optimum_hyperparameter.txt, SVR_optimum_hyperparameter.txt, and RF_optimum_hyperparameter.txt.
  2. We then tested the trained models in testing populations (MESA and METS). Now, because the non-linear ML models do not generate weight/beta values that we could have stored as DB files and used to predict expression in a new test set, the non-linear ML trained models in this context practically refers to the optimum hyperparameter files where we have stored the optimal hyperparameter values as referenced in step 1 above. Therefore, in order to use the non-linear ML models to predict expression in a new test set, firstly, for each gene, we essentially fit each of the non-linear ML algorithm with the training data as well as the optimal hyperparameter value learned for that gene in step 1 above, and then lastly we use the fitted model to predict expression in a new test set. This second step is actualized with 01_model_testing_in_METS.py, and 02_model_imputation_in_MESA.py scripts.
  
  # Model Usage and Caveats

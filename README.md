# ML-PredictDB
In this study we sought to optimize protein coding gene expression imputation performance across global populations in comparison to the popular transcriptome tool [PrediXcan](https://github.com/WheelerLab/PrediXcan). Specifically, we used non-linear machine learning (ML) models viz; random forest (RF), support vector regression (SVR), and K nearest neighbor (KNN) to build transcriptome imputation models, and evaluated their performance in comparison to elastic net (EN) - the algorithm used in [PrediXcan](https://github.com/WheelerLab/PrediXcan). We trained gene expression prediction models using genotype and blood monocyte transcriptome data from the Multi-Ethnic Study of Atherosclerosis (MESA) comprising individuals of African (AFA), Hispanic (HIS), and European (CAU) ancestries and tested them using genotype and whole blood transcriptome data from the Modeling the Epidemiology Transition Study (METS) comprising individuals of African ancestries. We show that the prediction performance is highest when the training and the testing population share similar ancestries regardless of the prediction algorithm used. While EN generally outperformed RF, SVR, and KNN, we found that RF outperforms EN for some genes, particularly between disparate ancestries, suggesting potential robustness and reduced variability of RF imputation performance across global populations. When applied to a high-density lipoprotein (HDL) phenotype, we show including RF prediction models in PrediXcan reveals potential gene associations missed by EN models. Therefore, by integrating non-linear modeling into PrediXcan and diversifying our training populations to include more global ancestries, we may uncover new genes associated with complex traits.

# Model Training and Testing
The model training and testing were done in two separate parts. Firstly, the elastic net (EN) models were trained and tested as described in [Mogil et al., 2018](https://github.com/WheelerLab/DivPop).
The non-linear machine learning (ML) models were built and tested as described below:
  1. We first used the 00_gridsearch.py script to carry out a gridsearch through each ML algorithm hyperparameters in order to determine the hyperparameter values or combination of values that yields the optimum five-fold cross-validated imputation performance (R<sup>2</sup>) for each gene. The optimum hyperparameter values and (R<sup>2</sup>) for each gene across training population and ML algorithm.

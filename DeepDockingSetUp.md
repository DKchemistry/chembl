Deep Docking Set Up 

(1) Pareto ranking to replace to docking score for enrichment metrics 

(2) Cheminformatic Analysis of likely contributors to strain 

(3) (LIT-PCBA)

(3.5) Deep Docking Set Up (What kind/how many inactives?)

(4) Optuna/Bayseian Optimization of Hyperparameters 

(5) Ensemble Model Approach, Aquisition Function

(6) Infastructure related upgrades TF1 - TF2 or going to a later PyTorch 

- cLog P vs MW to compare actives/inactive (like Irwin)

- Pareto in place of scoring functions 

- Molpal: Monte Carlo Drop out
- - Time consuming
- - Used to estimate uncertainity 
- - Not likely relevant to us (or even them) because we ultimately find "greedy" aquisition works the best (probability of best scoring)

- Our "greedy" approach is not super strict. 

- Hyperparameter Optimization via Optuna (!) 

- Multi-layer perceptron is likely to stay 

- Ensemble of models approach from Optuna (like concensus scoring in a sense) 

- However an ensemable of models is very similar Monte Carlo dropout 

- Our uncertainty from the classifier model is not the same as "ensemble" uncertainty, that is "epistemic uncertainty" (Francesco isn't sure about this)

- Maybe PyTorch model in place? Or Tensorflow 2? 
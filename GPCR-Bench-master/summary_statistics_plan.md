Need to get summary statistics
Need to get pareto statistics 

what summary statistics do we want? 

What questions do we want to answer? 

- What strain cut off yielded the best improvement in LogAUC?
- What were the baseline LogAUCs? 

`protein, baseline log auc, best strain log auc improvement` as a csv 

- what are the summary statistics within a protein for strain? 

`strain, ef1%, ef5%, logAUC`
 4, 3.5, 5.5, 0.24
4.5 ...

and so on for all the proteins 

we would need to get the ef1% and ef%5 per strain as well. 

then we could have some file like 

protein, ef1%, ef5%, logauc, strain threshold, strain ef1%, strain ef%5, strain logauc, pareto, pareto ef1%, pareto ef5%, pareto logauc 

maybe have pareto seperate? and baseline seperate? just have the strain as a multiline csv 

and really i have all of this data except ef1% and ef5% for strain, and pareto 

i think we just build this pipeline for like GPR40 and extend it all the way. 

we can keep the the same sort of input structure of files, and make that be case for lit pcba if we do that 

so the pipeline is 

input -> merged/concat dataframe with activity labels -> baseline performance (ef 1%, ef 5%, roc auc, semilogx auc) [output csv], pareto performance (ef 1%, ef 5%, roc auc, semilogx auc) [output csv], strain performance on range in thresholds (ef 1%, ef 5%, roc auc, semilogx auc) with delta performance in each metrics 

then those csvs can be combined, respectively, per protein. 3 csvs at the end (baseline, pareto, strain) (all proteins). those csvs can be processed all together for summary statistics

we can use the current GPR40.ipynb and start to clean it up for this task. i can show our work to avoid bad time management. 





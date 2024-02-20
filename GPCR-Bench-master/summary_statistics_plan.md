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

we can use the current GPR40.ipynb and start to clean it up for this task. i can show our work to avoid bad time management. we do not need to worry about the graphs for now, we can rework those seperately, but it is nice to have them as a sanity check. however the specific formatting is not very important 

okay i have the strain metrics for 1 protein, need baseline, and pareto.

let's do baseline and then think about other stuff 

# papermill

we need parameters for: 

title suffix

file processing (paths)

inactives or decoy

maybe style but not that important 

we should be able to accomodate multiple datasets - like gpcr bench and lit pcba, so maybe a title prefix 

maybe the pareto ranks themselves 

ultimately in psuedo code we want: 

make a notebook for every file proessing option as the name of title_prefix+title_suffix.ipynb 

it would be good to make the lit pcba dataset organized like gpcr bench is, but it seems kinda tricky to do that as they are structured differently in terms of data, so i think it might be best to have that as a fully seperate pipeline.


# pareto 

20 ranks is sufficient for francesco 

we need to write metrics for 

protein, Strain Energy Cutoff, EF1%, EF5%, ROC
xxx, Pareto 20,

Protein,Strain Energy Cutoff,EF1%,EF5%,deltaEF1%,deltaEF5%,Linear Log10 AUC (x10),Delta Linear Log10 AUC (x10),ROC_AUC,Actives,Total Count,deltaAUC


# Searching for Inactives in ChEMBL 

Using a well studied GPCR, Adenosine A2a receptor, [target report card.](https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL251/)

Trying ['Activities'](https://www.ebi.ac.uk/chembl/web_components/explore/activities/STATE_ID:1nnTnqAEceY1RZb6N9_drQ==)

How do I get the comment field? Previously I read this ['Using ChEMBL activity comments'](https://chembl.blogspot.com/2020/08/using-chembl-activity-comments.html) 

That link should be read carefully, the takeaways:  

1. Activity comments can refer to multiple things depending on how your search is done. 'Inactive' in an ADME assay is not the same as 'Inactive' in a functional assay. 
2. There are different data sources in ChEMBL, 'standard data' is what comes from the 7 core journals. I would *guess* that an IC50 of 0 for some 'standard' sourced data does mean 0, but I am not sure. 

In the examples from the link, note that the 'Standard Type' field is changing. A confusing (but likely relevant case) is the second/third image, where both the words used in the standard type (Activity, Potency) are then assigned comments that mean the same thing, but use different words (Not active, inactive)

This makes me confused about how to best search for inactive compounds for some GPCR - how many possible standard types are there and how many possible variations of the word inactive are there? I am guessing that is why the LIT PCBA paper went with PubChem data. Using their workflow for inactive selection would be ideal, but I don't see any of their code if they queried API's programmatically. Here is what they say. 

> Data Set
>Bioactivity data were retrieved from the PubChem BioAssay database, (22) where all information on true active and true inactive substances for a protein target is provided based on experimental results from confirmatory dose–response bioactivity assays, whose related details including assay principles, general protocols, and other remarks are also given. All data were updated as of December 31, 2018. The “limits” search engine (https://www.ncbi.nlm.nih.gov/pcassay/limits) was used to filter the PubChem BioAssays resource by various options, with “Activity Outcome” set as “Active”, “Substance Type” set as “Chemical”, and “Screening Stage” defined as “Confirmatory, Dose-Response”. Here, 149 assays targeting a single protein target, operated on at least 10,000 substances, and leading to at least 50 confirmed actives were first retained. The experimental screening data were kept if the target was characterized by at least one Protein Data Bank (PDB) (34) entry, in complex with a ligand of the same phenotype (i.e., inhibitor, agonist, or antagonist) as that of the tested active substances of the corresponding bioactivity assay. Altogether, 21 raw HTS data tables were directly retrieved as csv files from the PubChem BioAssay website as well as actives and inactives in separate sd files. The PDB resource was then browsed by Uniprot identifier (Uniprot ID) (35) to retrieve the corresponding PDB entries in the suitable ligand-bound form

So one strategy then is to try this approach and see if something like Beta2 antagonists are there. I've tried repeating their method as described here but getting null results and errors. 

Thomas wanted more actives, so for the sake of having actionable results, I am going to switch gears to searching for actives.
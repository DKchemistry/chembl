# Searching for Actives: OPRK1

## GPCR Bench Data

Let's start by finding actives for a dataset where we already have known inactives like OPRK1. The other benefit is that there been many newly published OPRK1 structures.  

Then, we can simply add our newly found actives, deduplicate them, and likely dock them without even having access to CCU clusters.  

Let's use the criteria from GPCR Bench first, they mentioned that they plan to automatically update GPCR bench. I checked. The files appear to be on the DUD-E website but the downloads either fail or are corrupt.  Found their github and tried there. The [github is working](https://github.com/dahliaweiss/GPCR-Bench/blob/master/OPRK_chembl19_set.smi) but it is only the the old data from 9 years ago. Maybe worthwhile to email the authors about either the code, but for now - I can at least get the actives data.  

I am getting the zip from github, in the meanwhile - the *.smi files contain the whole set post generation. The names are also given in an .id file.

```sh
tail -n 201 OPRK_chembl19_set.smi > bench_actives_oprk1.smi
```

The 201 ligand is the cocrystal ligand.

We need to check overlap with LIT PCBA OPRK1. I moved the .smi file to here as 'litpcba_actives_oprk1.smi'.

Will ICM be able to check identical mols? There might an issue with charge states, as GPCR Bench allowed for multiple protomers and their smiles are enumerated.

```sh
icm/def> find table exact litpcba_actives_oprk1 bench_actives_oprk1.mol select 0.5 ( no ? : charge ) 
 Info> 0 hits selected in litpcba_actives_oprk1
icm/def> find table exact bench_actives_oprk1 litpcba_actives_oprk1.mol select 0.5 ( no ? : charge ) 
 Info> 0 hits selected in bench_actives_oprk1
```

comparison within the set also looks fine

```sh
icm/def> find table bench_actives_oprk1 select index = Index( bench_actives_oprk1.mol exact ) 
 Info> 0 hits selected in bench_actives_oprk1
```

```sh
icm/def> find table litpcba_actives_oprk1 select index = Index( litpcba_actives_oprk1.mol exact ) 
 Info> 0 hits selected in litpcba_actives_oprk1
 ```

 Okay, so, at the least we have another 200 actives.  

 Since we should already have the data for lit pcba, we just have to prepare these and dock them.

 Unfortunately, I think my receptor is not on this laptop?  

 Either way:

 1. Rdkit
 2. Ligprep with no additional protonation
 3. Dock into grid
 4. Strain
 5. merge dfs
 6. performance metrics

## Initial ChEMBL Workflow

I think it is best to just get my automation scripts back and just get actives for the receptors thomas wanted.

We can return back to OPRK1 and work towards others.  

Let's a model a workflow, an initial search of `OPRK1` returns this [link](https://www.ebi.ac.uk/chembl/g/#search_results/targets/query=OPRK1)

There are already a lot of options to chose from. There are choices about species and 'groups' (I would guess this is for non selective molecules). I am not sure exactly how we would want to filter these. Protein species should probably be tied to our PDB availability, but practically - we want to use human PDBs. Selectivity is a bit more nebulous but I would avoid stuff like this:

| ChEMBL ID | Name | UniProt Accessions | Type |
| --- | --- | --- | --- |
| CHEMBL2095192 | Opioid receptor | P42866\|P32300\|P33534 | PROTEIN FAMILY |

Instead we'll look at this:

| ChEMBL ID | Name | UniProt Accessions | Type | Organism | Compounds |
| --- | --- | --- | --- | --- | --- |
| [CHEMBL237](https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL237/) | Kappa opioid receptor | P41145 | SINGLE PROTEIN | Homo sapiens | 7926 |

So now, we'll follow the [CHEMBL237 link](https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL237/). This is the target card information. It takes a very long time to load.  

However, all we really want is probably the [bioactivities](https://www.ebi.ac.uk/chembl/web_components/explore/activities/STATE_ID:D3aRKqliHl08KJtnnD4VfA==). We get 14,710 records. How do we filter this down? These are all the options available.

molecule_chembl_id | Molecule ChEMBL ID | string |
| --- | --- | --- |
| molecule_pref_name | Molecule Name | string |
| _metadata.parent_molecule_data.max_phase | Molecule Max Phase | string |
| _metadata.parent_molecule_data.full_mwt | Molecular Weight | double |
| _metadata.parent_molecule_data.num_ro5_violations | #RO5 Violations | integer |
| _metadata.parent_molecule_data.alogp | AlogP | double |
| _metadata.parent_molecule_data.compound_key | Compound Key | string |
| canonical_smiles | Smiles | string |
| standard_type | Standard Type | string |
| standard_relation | Standard Relation | string |
| standard_value | Standard Value | double |
| standard_units | Standard Units | string |
| pchembl_value | pChEMBL Value | double |
| data_validity_comment | Data Validity Comment | string |
| activity_comment | Comment | string |
| uo_units | Uo Units | string |
| ligand_efficiency.bei | Ligand Efficiency BEI | double |
| ligand_efficiency.le | Ligand Efficiency LE | double |
| ligand_efficiency.lle | Ligand Efficiency LLE | double |
| ligand_efficiency.sei | Ligand Efficiency SEI | double |
| potential_duplicate | Potential Duplicate | integer |
| assay_chembl_id | Assay ChEMBL ID | string |
| assay_description | Assay Description | string |
| assay_type | Assay Type | string |
| bao_format | BAO Format ID | string |
| bao_label | BAO Label | string |
| _metadata.assay_data.assay_organism | Assay Organism | string |
| _metadata.assay_data.tissue_chembl_id | Assay Tissue ChEMBL ID | string |
| _metadata.assay_data.assay_tissue | Assay Tissue Name | string |
| _metadata.assay_data.assay_cell_type | Assay Cell Type | string |
| _metadata.assay_data.assay_subcellular_fraction | Assay Subcellular Fraction | string |
| _metadata.assay_data.assay_parameters | Assay Parameters | object |
| assay_variant_accession | Assay Variant Accession | string |
| assay_variant_mutation | Assay Variant Mutation | string |
| target_chembl_id | Target ChEMBL ID | string |
| target_pref_name | Target Name | string |
| target_organism | Target Organism | string |
| _metadata.target_data.target_type | Target Type | string |
| document_chembl_id | Document ChEMBL ID | string |
| src_id | Source ID | string |
| _metadata.source.src_description | Source Description | string |
| document_journal | Document Journal | string |
| document_year | Document Year | integer |
| _metadata.assay_data.cell_chembl_id | Cell ChEMBL ID | string |
| activity_properties | Properties | object |
| action_type.action_type | Action Type | string |

And for reference, this is the query if we do need to move to something programmatic: 

Index Name  
`chembl_activity`

```json
{
  "query": {
    "bool": {
      "must": [
        {
          "query_string": {
            "query": "target_chembl_id:CHEMBL237"
          }
        }
      ],
      "filter": [
        {
          "query_string": {
            "analyze_wildcard": true,
            "query": "target_chembl_id:CHEMBL237"
          }
        }
      ],
      "should": [],
      "must_not": []
    }
  },
  "track_total_hits": true,
  "size": 20,
  "from": 0,
  "_source": [
    "molecule_chembl_id",
    "_metadata.parent_molecule_data.compound_key",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "activity_comment",
    "assay_chembl_id",
    "assay_description",
    "bao_label",
    "_metadata.assay_data.assay_organism",
    "target_chembl_id",
    "target_pref_name",
    "target_organism",
    "_metadata.target_data.target_type",
    "document_chembl_id",
    "_metadata.source.src_description",
    "_metadata.assay_data.cell_chembl_id",
    "molecule_pref_name",
    "_metadata.parent_molecule_data.max_phase",
    "_metadata.parent_molecule_data.full_mwt",
    "_metadata.parent_molecule_data.num_ro5_violations",
    "_metadata.parent_molecule_data.alogp",
    "canonical_smiles",
    "data_validity_comment",
    "uo_units",
    "ligand_efficiency.bei",
    "ligand_efficiency.le",
    "ligand_efficiency.lle",
    "ligand_efficiency.sei",
    "potential_duplicate",
    "assay_type",
    "bao_format",
    "_metadata.assay_data.tissue_chembl_id",
    "_metadata.assay_data.assay_tissue",
    "_metadata.assay_data.assay_cell_type",
    "_metadata.assay_data.assay_subcellular_fraction",
    "_metadata.assay_data.assay_parameters",
    "assay_variant_accession",
    "assay_variant_mutation",
    "src_id",
    "document_journal",
    "document_year",
    "activity_properties",
    "action_type.action_type"
  ],
  "sort": []
}
```

Here are there filters:

| Property                    | Filter                  |
|-----------------------------|-------------------------|
| Ki, Kd, IC50, or EC50       | <= 10,000 nM   |
| Molecular Weight (MW)       | 250 ≤ MW ≤ 500          |
| LogP                        | -2 ≤ LogP ≤ 5           |
| Hydrogen bond donors (HBD)  | HBD ≤ 4                 |
| Hydrogen bond acceptors (HBA)| 8 ≤ HBA ≤ 10            |
| Rotatable bonds             | Rotatable bonds < 8     |
| Rings                       | Rings < 4               |
| 9-membered rings            | 0                       |
| Isotopes                    | 0                       |
| Sulfur (S) atoms            | S ≤ 3                   |
| Chlorine (Cl) atoms         | Cl ≤ 3                  |
| Bromine (Br) atoms          | Br ≤ 2                  |
| Iodine (I) atoms            | I ≤ 2                   |
| Fluorine (F) atoms          | F ≤ 7                   |
| Heteroatoms                 | Only S, F, Cl, Br, I, N, O |
| MW < 250                    | LogP < 4.5              |
| PAINS                       | Filter Pains    |
| Functional groups           | "Specific filters apply"  |

Not all of this is possible to filter in ChEMBL as far as I can tell. There should be a way to do this relatively easily in ICM. I would probably just try to filter the activity data as such before hand. 

After 


| Property                    | Filter                  |
|-----------------------------|-------------------------|
| Ki, Kd, IC50, or EC50       | <= 10,000 nM   |
| LogP                        | -2 ≤ LogP ≤ 5           |

We are already down to 5,452. I grabbed the TSV/CSV now as the site seems unstable. Doing this in ICM without dark mode would burn my eyes out of my skull. I can do that in the morning. The other thing we should filter are the two data sets they say are unreliable.

Here is the query just in case:

```json
{
  "query": {
    "bool": {
      "must": [
        {
          "query_string": {
            "query": "target_chembl_id:CHEMBL237"
          }
        }
      ],
      "filter": [
        {
          "query_string": {
            "analyze_wildcard": true,
            "query": "target_chembl_id:CHEMBL237"
          }
        },
        {
          "bool": {
            "should": [
              {
                "terms": {
                  "standard_type": [
                    "Ki",
                    "EC50",
                    "IC50"
                  ]
                }
              }
            ]
          }
        },
        {
          "range": {
            "_metadata.parent_molecule_data.alogp": {
              "gte": -2,
              "lt": 5
            }
          }
        },
        {
          "range": {
            "standard_value": {
              "gte": -56,
              "lt": "10000"
            }
          }
        }
      ],
      "should": [],
      "must_not": []
    }
  },
  "track_total_hits": true,
  "size": 20,
  "from": 0,
  "_source": [
    "molecule_chembl_id",
    "_metadata.parent_molecule_data.compound_key",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "activity_comment",
    "assay_chembl_id",
    "assay_description",
    "bao_label",
    "_metadata.assay_data.assay_organism",
    "target_chembl_id",
    "target_pref_name",
    "target_organism",
    "_metadata.target_data.target_type",
    "document_chembl_id",
    "_metadata.source.src_description",
    "_metadata.assay_data.cell_chembl_id",
    "molecule_pref_name",
    "_metadata.parent_molecule_data.max_phase",
    "_metadata.parent_molecule_data.full_mwt",
    "_metadata.parent_molecule_data.num_ro5_violations",
    "_metadata.parent_molecule_data.alogp",
    "canonical_smiles",
    "data_validity_comment",
    "uo_units",
    "ligand_efficiency.bei",
    "ligand_efficiency.le",
    "ligand_efficiency.lle",
    "ligand_efficiency.sei",
    "potential_duplicate",
    "assay_type",
    "bao_format",
    "_metadata.assay_data.tissue_chembl_id",
    "_metadata.assay_data.assay_tissue",
    "_metadata.assay_data.assay_cell_type",
    "_metadata.assay_data.assay_subcellular_fraction",
    "_metadata.assay_data.assay_parameters",
    "assay_variant_accession",
    "assay_variant_mutation",
    "src_id",
    "document_journal",
    "document_year",
    "activity_properties",
    "action_type.action_type"
  ],
  "sort": []
}
```

I can apply this same sort of logic to whatever GPCRs we want. As long as Thomas agrees we should be fine. In the morning we can discuss 

1. Decoys or known inactives? Both?
2. What ratio of actives/inactives? 
3. What makes a good benchmark? Both a "worst case" and a "best case" seem reasonable. A prospective study *definitely* warrants a best case. 
  
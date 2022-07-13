#!/bin/bash
(
cat <<HERE
{
  "filters":{
    "op":"and",
    "content":[
      {
        "op":"in",
        "content":{
          "field":"file_name",
          "value":[
HERE
tail -n +2 "$1" | awk '{print $2}' | perl -pe 's/(.*)/"\1"/' | paste -s -d ','
cat <<HERE
          ]
        }
      }
    ]
  },
  "format":"tsv",
  "fields":"file_id,file_name,experimental_strategy,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id,cases.demographic.race,cases.demographic.ethnicity,cases.demographic.gender",
  "size":"40000"
}
HERE
) | curl --request POST --header 'Content-Type: application/json' --data '@-' 'https://api.gdc.cancer.gov/files'

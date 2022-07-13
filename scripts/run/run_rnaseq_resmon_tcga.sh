#!/bin/bash
# Start carefully. We can adapt the job number dynamically.
echo 1 > /tmp/rna_jobs_resmon_tcga

parallel --progress \
--limit 'scripts/tools/scripts/resource_limit.sh' \
-j /tmp/rna_jobs_resmon_tcga \
bash scripts/tools/scripts/adapter_tcga_sample.sh '{2}' '{1}' \
'temp/output_resmon_tcga' \
'temp/resmon_single_tcga3' \
'temp/nosnapshot/log_tcga_resmon_rna' \
'$PARALLEL_JOBSLOT' :::: \
'temp/resmon_single_tcga3_rna_files.tsv' \
<(ls scripts/tools/RNASeq/*.sh )

echo "Finished"
#!/bin/bash
# Start carefully. We can adapt the job number dynamically.
echo 4 > /tmp/rna_jobs_tcga

parallel --progress \
--limit 'scripts/tools/scripts/resource_limit.sh' \
-j /tmp/rna_jobs_tcga \
bash scripts/tools/scripts/adapter_tcga_sample.sh '{2}' '{1}' \
'temp/output_tools_tcga' \
'temp/resmon_tcga_dummy' \
'temp/nosnapshot/log_tcga_rna' \
'$PARALLEL_JOBSLOT' :::: \
'temp/tcga_bam_list_rna.tsv' \
<(ls scripts/tools/RNASeq/*.sh )

echo "Finished"
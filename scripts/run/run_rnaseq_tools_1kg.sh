#!/bin/bash
# Start carefully. We can adapt the job number dynamically.
echo 4 > /tmp/rna_jobs_tcga

parallel --progress \
--limit 'scripts/tools/scripts/resource_limit.sh' \
-j /tmp/rna_jobs_tcga \
bash scripts/tools/scripts/adapter.sh '{2}' '{1}' \
'temp/output_tools_1kg' \
'temp/resmon_1kg_dummy' \
'temp/nosnapshot/log_1kg_rna' \
'$PARALLEL_JOBSLOT' :::: \
<(find 'downloads/1kg/bam_geuvadis/' -maxdepth 1 -name '*.bam' ) \
<(ls scripts/tools/RNASeq/*.sh )

echo "Finished"
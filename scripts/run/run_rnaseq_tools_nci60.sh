#!/bin/bash
# Start carefully. We can adapt the job number dynamically.
echo 4 > /tmp/rna_jobs_nci60

parallel --progress \
--limit 'scripts/tools/scripts/resource_limit.sh' \
-j /tmp/rna_jobs_nci60 \
bash scripts/tools/scripts/adapter.sh '{2}' '{1}' \
'temp/output_tools_nci60_1kg' \
'temp/resmon_nci60_dummy' \
'temp/nosnapshot/logrna_nci60' \
'$PARALLEL_JOBSLOT' :::: \
<(find downloads/pub/reinhold_2019/HLA/ -name '*.bam' ) \
<(ls scripts/tools/RNASeq/*.sh )

echo "Finished"
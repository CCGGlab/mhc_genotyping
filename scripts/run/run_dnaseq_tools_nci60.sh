#!/bin/bash
# Start carefully. We can adapt the job number dynamically.
echo 4 > /tmp/dna_jobs_nci60

parallel --progress \
--limit 'scripts/tools/scripts/resource_limit.sh' \
-j /tmp/dna_jobs_nci60 \
bash scripts/tools/scripts/adapter.sh '{2}' '{1}' \
'temp/output_tools_nci60_1kg' \
'temp/resmon_nci60_dummy' \
'temp/nosnapshot/logwes_nci60' \
'$PARALLEL_JOBSLOT' :::: \
<(find 'downloads/pub/abaan_2013/HLA' -name '*.bam' ) \
<(ls scripts/tools/DNASeq/*.sh )

echo "Finished"
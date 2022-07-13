#!/bin/bash
# Start carefully. We can adapt the job number dynamically.
echo 4 > /tmp/dna_jobs_1kg

parallel --progress \
--limit 'scripts/tools/scripts/resource_limit.sh' \
-j /tmp/dna_jobs_1kg \
bash scripts/tools/scripts/adapter.sh '{2}' '{1}' \
'temp/output_tools_1kg' \
'temp/resmon_1kg_dummy' \
'temp/nosnapshot/log_1kg_dna' \
'$PARALLEL_JOBSLOT' :::: \
<(find 'downloads/1kg/bam_hs38/merged_alt_mapped/' -maxdepth 1 -name '*.bam' ) \
<(ls scripts/tools/DNASeq/*.sh )

echo "Finished"
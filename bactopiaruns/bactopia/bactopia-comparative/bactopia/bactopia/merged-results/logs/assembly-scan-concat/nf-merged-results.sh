#!/bin/bash -ue
csvtk \
    concat \
     \
    --num-cpus 2 \
    --tabs  \
    --out-tabs \
    --out-file assembly-scan.tsv \
    SRR2838702.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ASSEMBLER:CSVTK_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

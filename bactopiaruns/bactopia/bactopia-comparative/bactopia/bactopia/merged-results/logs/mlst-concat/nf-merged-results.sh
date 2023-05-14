#!/bin/bash -ue
csvtk \
    concat \
    --no-header-row \
    --num-cpus 2 \
    --tabs  \
    --out-tabs \
    --out-file mlst.tsv \
    SRR2838702.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:MLST:CSVTK_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

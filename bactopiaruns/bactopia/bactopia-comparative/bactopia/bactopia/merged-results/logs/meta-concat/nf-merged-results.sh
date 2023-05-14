#!/bin/bash -ue
csvtk \
    concat \
     \
    --num-cpus 2 \
    --tabs  \
    --out-tabs \
    --out-file meta.tsv \
    SRR2838702-meta.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:GATHER:CSVTK_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

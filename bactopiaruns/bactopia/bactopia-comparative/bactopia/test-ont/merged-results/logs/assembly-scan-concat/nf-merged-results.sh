#!/bin/bash -ue
csvtk \
    concat \
     \
    --num-cpus 4 \
    --tabs  \
    --out-tabs \
    --out-file assembly-scan.tsv \
    test-ont.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ASSEMBLER:CSVTK_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

#!/bin/bash -ue
csvtk \
    concat \
    --no-header-row \
    --num-cpus 4 \
    --tabs  \
    --out-tabs \
    --out-file mlst.tsv \
    test-pe.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:MLST:CSVTK_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

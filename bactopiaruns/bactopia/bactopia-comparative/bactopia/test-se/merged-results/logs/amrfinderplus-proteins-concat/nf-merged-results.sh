#!/bin/bash -ue
csvtk \
    concat \
     \
    --num-cpus 4 \
    --tabs  \
    --out-tabs \
    --out-file amrfinderplus-proteins.tsv \
    test-se-proteins.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:PROTEINS_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

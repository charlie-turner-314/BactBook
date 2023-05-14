#!/bin/bash -ue
csvtk \
    concat \
     \
    --num-cpus 4 \
    --tabs  \
    --out-tabs \
    --out-file amrfinderplus-genes.tsv \
    test-se-genes.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:GENES_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

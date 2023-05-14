#!/bin/bash -ue
csvtk \
    concat \
     \
    --num-cpus 2 \
    --tabs  \
    --out-tabs \
    --out-file amrfinderplus-genes.tsv \
    SRR2838702-genes.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:GENES_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

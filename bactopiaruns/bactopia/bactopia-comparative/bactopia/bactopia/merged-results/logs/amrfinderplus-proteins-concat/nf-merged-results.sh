#!/bin/bash -ue
csvtk \
    concat \
     \
    --num-cpus 2 \
    --tabs  \
    --out-tabs \
    --out-file amrfinderplus-proteins.tsv \
    SRR2838702-proteins.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:PROTEINS_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS

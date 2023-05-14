#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-se.ffn.gz > test-se.ffn
fi

if [ "true" == "true" ]; then
    gzip -c -d test-se.faa.gz > test-se.faa
fi

tar xzvf amrfinderplus.tar.gz

# Gene
amrfinder \
   -n test-se.ffn \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-se > test-se-genes.tsv

# Protein
amrfinder \
   -p test-se.faa \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-se > test-se-proteins.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:AMRFINDERPLUS_RUN":
    amrfinderplus: $(amrfinder --version)
    amrfinderplus-database: $(echo $(echo $(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
END_VERSIONS

#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-ont.ffn.gz > test-ont.ffn
fi

if [ "true" == "true" ]; then
    gzip -c -d test-ont.faa.gz > test-ont.faa
fi

tar xzvf amrfinderplus.tar.gz

# Gene
amrfinder \
   -n test-ont.ffn \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-ont > test-ont-genes.tsv

# Protein
amrfinder \
   -p test-ont.faa \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-ont > test-ont-proteins.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:AMRFINDERPLUS_RUN":
    amrfinderplus: $(amrfinder --version)
    amrfinderplus-database: $(echo $(echo $(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
END_VERSIONS

#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-pe.ffn.gz > test-pe.ffn
fi

if [ "true" == "true" ]; then
    gzip -c -d test-pe.faa.gz > test-pe.faa
fi

tar xzvf amrfinderplus.tar.gz

# Gene
amrfinder \
   -n test-pe.ffn \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-pe > test-pe-genes.tsv

# Protein
amrfinder \
   -p test-pe.faa \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-pe > test-pe-proteins.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:AMRFINDERPLUS_RUN":
    amrfinderplus: $(amrfinder --version)
    amrfinderplus-database: $(echo $(echo $(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
END_VERSIONS

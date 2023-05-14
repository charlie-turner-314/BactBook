#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-short-polish.ffn.gz > test-short-polish.ffn
fi

if [ "true" == "true" ]; then
    gzip -c -d test-short-polish.faa.gz > test-short-polish.faa
fi

tar xzvf amrfinderplus.tar.gz

# Gene
amrfinder \
   -n test-short-polish.ffn \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-short-polish > test-short-polish-genes.tsv

# Protein
amrfinder \
   -p test-short-polish.faa \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 4 \
    --name test-short-polish > test-short-polish-proteins.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:AMRFINDERPLUS_RUN":
    amrfinderplus: $(amrfinder --version)
    amrfinderplus-database: $(echo $(echo $(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
END_VERSIONS

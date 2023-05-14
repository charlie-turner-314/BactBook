#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d SRR2838702.ffn.gz > SRR2838702.ffn
fi

if [ "true" == "true" ]; then
    gzip -c -d SRR2838702.faa.gz > SRR2838702.faa
fi

tar xzvf amrfinderplus.tar.gz

# Gene
amrfinder \
   -n SRR2838702.ffn \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 2 \
    --name SRR2838702 > SRR2838702-genes.tsv

# Protein
amrfinder \
   -p SRR2838702.faa \
     \
    --plus --ident_min -1 --coverage_min 0.5 --translation_table 11 \
    --database amrfinderplus/ \
    --threads 2 \
    --name SRR2838702 > SRR2838702-proteins.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:AMRFINDERPLUS:AMRFINDERPLUS_RUN":
    amrfinderplus: $(amrfinder --version)
    amrfinderplus-database: $(echo $(echo $(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
END_VERSIONS

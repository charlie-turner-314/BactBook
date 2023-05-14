#!/bin/bash -ue
tar -xzvf mlst.tar.gz

mlst \
    --threads 4 \
    --blastdb mlstdb/blast/mlst.fa \
    --datadir mlstdb/pubmlst \
    --minid 95 --mincov 10 --minscore 50 \
    test-pe.fna.gz \
    > test-pe.tsv

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:MLST:MLST_MODULE":
    mlst: $( echo $(mlst --version 2>&1) | sed 's/mlst //' )
    mlst-database: $( cat mlstdb/DB_VERSION )
END_VERSIONS

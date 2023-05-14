#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-pe.fna.gz > test-pe.fna
fi

prokka \
    --evalue 1e-09 --coverage 80 --centre Bactopia \
    --cpus 4 \
    --prefix test-pe \
     \
    --locustag test-pe \
    --proteins proteins.faa \
     \
    test-pe.fna

# Make blastdb of contigs, genes, proteins
mkdir blastdb
cat test-pe/test-pe.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for !{meta.id}" -out blastdb/!{meta.id}.fna
cat test-pe/test-pe.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for !{meta.id}" -out blastdb/!{meta.id}.ffn
cat test-pe/test-pe.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for !{meta.id}" -out blastdb/!{meta.id}.faa
tar -cvf - blastdb/ | gzip -c > test-pe/test-pe-blastdb.tar.gz

if [[ "false" == "false" ]]; then
    gzip test-pe/*.gff
    gzip test-pe/*.gbk
    gzip test-pe/*.fna
    gzip test-pe/*.faa
    gzip test-pe/*.ffn
    gzip test-pe/*.sqn
    gzip test-pe/*.fsa
    gzip test-pe/*.tbl
fi

mv test-pe/test-pe.err ./
mv test-pe/test-pe.log ./
mv test-pe/ results/

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ANNOTATOR:PROKKA_MODULE":
    makeblastdb: $( echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    prokka: $( echo $(prokka --version 2>&1) | sed 's/^.*prokka //')
END_VERSIONS

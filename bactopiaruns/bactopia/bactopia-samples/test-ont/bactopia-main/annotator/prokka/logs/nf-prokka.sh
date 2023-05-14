#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-ont.fna.gz > test-ont.fna
fi

prokka \
    --evalue 1e-09 --coverage 80 --centre Bactopia \
    --cpus 4 \
    --prefix test-ont \
     \
    --locustag test-ont \
    --proteins proteins.faa \
     \
    test-ont.fna

# Make blastdb of contigs, genes, proteins
mkdir blastdb
cat test-ont/test-ont.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for !{meta.id}" -out blastdb/!{meta.id}.fna
cat test-ont/test-ont.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for !{meta.id}" -out blastdb/!{meta.id}.ffn
cat test-ont/test-ont.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for !{meta.id}" -out blastdb/!{meta.id}.faa
tar -cvf - blastdb/ | gzip -c > test-ont/test-ont-blastdb.tar.gz

if [[ "false" == "false" ]]; then
    gzip test-ont/*.gff
    gzip test-ont/*.gbk
    gzip test-ont/*.fna
    gzip test-ont/*.faa
    gzip test-ont/*.ffn
    gzip test-ont/*.sqn
    gzip test-ont/*.fsa
    gzip test-ont/*.tbl
fi

mv test-ont/test-ont.err ./
mv test-ont/test-ont.log ./
mv test-ont/ results/

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ANNOTATOR:PROKKA_MODULE":
    makeblastdb: $( echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    prokka: $( echo $(prokka --version 2>&1) | sed 's/^.*prokka //')
END_VERSIONS

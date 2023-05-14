#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-se.fna.gz > test-se.fna
fi

prokka \
    --evalue 1e-09 --coverage 80 --centre Bactopia \
    --cpus 4 \
    --prefix test-se \
     \
    --locustag test-se \
    --proteins proteins.faa \
     \
    test-se.fna

# Make blastdb of contigs, genes, proteins
mkdir blastdb
cat test-se/test-se.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for !{meta.id}" -out blastdb/!{meta.id}.fna
cat test-se/test-se.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for !{meta.id}" -out blastdb/!{meta.id}.ffn
cat test-se/test-se.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for !{meta.id}" -out blastdb/!{meta.id}.faa
tar -cvf - blastdb/ | gzip -c > test-se/test-se-blastdb.tar.gz

if [[ "false" == "false" ]]; then
    gzip test-se/*.gff
    gzip test-se/*.gbk
    gzip test-se/*.fna
    gzip test-se/*.faa
    gzip test-se/*.ffn
    gzip test-se/*.sqn
    gzip test-se/*.fsa
    gzip test-se/*.tbl
fi

mv test-se/test-se.err ./
mv test-se/test-se.log ./
mv test-se/ results/

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ANNOTATOR:PROKKA_MODULE":
    makeblastdb: $( echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    prokka: $( echo $(prokka --version 2>&1) | sed 's/^.*prokka //')
END_VERSIONS

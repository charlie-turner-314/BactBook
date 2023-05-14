#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d test-short-polish.fna.gz > test-short-polish.fna
fi

prokka \
    --evalue 1e-09 --coverage 80 --centre Bactopia \
    --cpus 4 \
    --prefix test-short-polish \
     \
    --locustag test-short-polish \
    --proteins proteins.faa \
     \
    test-short-polish.fna

# Make blastdb of contigs, genes, proteins
mkdir blastdb
cat test-short-polish/test-short-polish.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for !{meta.id}" -out blastdb/!{meta.id}.fna
cat test-short-polish/test-short-polish.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for !{meta.id}" -out blastdb/!{meta.id}.ffn
cat test-short-polish/test-short-polish.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for !{meta.id}" -out blastdb/!{meta.id}.faa
tar -cvf - blastdb/ | gzip -c > test-short-polish/test-short-polish-blastdb.tar.gz

if [[ "false" == "false" ]]; then
    gzip test-short-polish/*.gff
    gzip test-short-polish/*.gbk
    gzip test-short-polish/*.fna
    gzip test-short-polish/*.faa
    gzip test-short-polish/*.ffn
    gzip test-short-polish/*.sqn
    gzip test-short-polish/*.fsa
    gzip test-short-polish/*.tbl
fi

mv test-short-polish/test-short-polish.err ./
mv test-short-polish/test-short-polish.log ./
mv test-short-polish/ results/

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ANNOTATOR:PROKKA_MODULE":
    makeblastdb: $( echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    prokka: $( echo $(prokka --version 2>&1) | sed 's/^.*prokka //')
END_VERSIONS

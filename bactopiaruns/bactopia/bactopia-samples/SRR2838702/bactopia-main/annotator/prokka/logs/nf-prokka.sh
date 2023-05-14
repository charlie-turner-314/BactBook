#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d SRR2838702.fna.gz > SRR2838702.fna
fi

prokka \
    --evalue 1e-09 --coverage 80 --centre Bactopia \
    --cpus 2 \
    --prefix SRR2838702 \
     \
    --locustag SRR2838702 \
    --proteins proteins.faa \
     \
    SRR2838702.fna

# Make blastdb of contigs, genes, proteins
mkdir blastdb
cat SRR2838702/SRR2838702.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for !{meta.id}" -out blastdb/!{meta.id}.fna
cat SRR2838702/SRR2838702.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for !{meta.id}" -out blastdb/!{meta.id}.ffn
cat SRR2838702/SRR2838702.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for !{meta.id}" -out blastdb/!{meta.id}.faa
tar -cvf - blastdb/ | gzip -c > SRR2838702/SRR2838702-blastdb.tar.gz

if [[ "false" == "false" ]]; then
    gzip SRR2838702/*.gff
    gzip SRR2838702/*.gbk
    gzip SRR2838702/*.fna
    gzip SRR2838702/*.faa
    gzip SRR2838702/*.ffn
    gzip SRR2838702/*.sqn
    gzip SRR2838702/*.fsa
    gzip SRR2838702/*.tbl
fi

mv SRR2838702/SRR2838702.err ./
mv SRR2838702/SRR2838702.log ./
mv SRR2838702/ results/

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ANNOTATOR:PROKKA_MODULE":
    makeblastdb: $( echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    prokka: $( echo $(prokka --version 2>&1) | sed 's/^.*prokka //')
END_VERSIONS

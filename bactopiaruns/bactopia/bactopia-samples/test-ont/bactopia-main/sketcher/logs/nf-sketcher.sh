#!/bin/bash -ue
if [ "true" == "true" ]; then
    xz -c -d mash-refseq88.k21.msh.xz > mash-refseq88.k21.msh
fi

gzip -cd test-ont.fna.gz | mash sketch -o test-ont-k21 -k 21 -s 10000 -I test-ont -
gzip -cd test-ont.fna.gz | mash sketch -o test-ont-k31 -k 31 -s 10000 -I test-ont -
sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=10000 --merge test-ont -o test-ont.sig test-ont.fna.gz

# Mash Screen

echo "identity<TAB>shared-hashes<TAB>median-multiplicity<TAB>p-value<TAB>query-ID<TAB>query-comment" | sed 's/<TAB>/	/g' > test-ont-mash-refseq88-k21.txt
gzip -cd test-ont.fna.gz | mash screen -w -i 0.8 -p 4 mash-refseq88.k21.msh.xz - | sort -gr >> test-ont-mash-refseq88-k21.txt

# Sourmash classify
sourmash lca classify --query test-ont.sig --db gtdb-rs207.genomic-reps.dna.k31.lca.json.gz > test-ont-sourmash-gtdb-rs207-k31.txt

# Capture versions
cat <<-END_VERSIONS > versions.yml
"BACTOPIA:SKETCHER:SKETCHER_MODULE":
    mash: $(echo $(mash 2>&1) | sed 's/^.*Mash version //;s/ .*$//')
    sourmash: $(echo $(sourmash --version 2>&1) | sed 's/sourmash //;')
END_VERSIONS

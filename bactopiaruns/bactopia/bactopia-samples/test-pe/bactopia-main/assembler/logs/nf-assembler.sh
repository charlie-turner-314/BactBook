#!/bin/bash -ue
if [ "null" == "true" ]; then
    mkdir results
    gzip -cd empty.fna.gz > results/test-pe.fna
elif [[ "paired-end" == "hybrid" || "false" == "true" ]]; then
    # Unicycler
    unicycler \
        -1 test-pe_R1.fastq.gz -2 test-pe_R2.fastq.gz \
        --keep 1 --min_fasta_length 500 --mode normal --min_component_size 1000 --min_dead_end_size 1000 \
         \
        -o results/ \
        --threads 4
    sed -r 's/^>([0-9]+)(.*)/>test-pe_\1\2/' results/assembly.fasta > results/test-pe.fna
    mv results/assembly.gfa results/unicycler.gfa
elif [[ "paired-end" == "ont" || "paired-end" == "short_polish" ]]; then
    # Dragonflye
    if ! dragonflye \
        --reads test-pe_R1.fastq.gz \
        --gsize 358000 \
        --outdir results \
        --assembler flye --polypolish 1 --minlen 500 --mincov 2 --force --keepfiles --depth 0 --minreadlen 0 --minquality 0 --racon 1 --medaka 0 \
        --namefmt "test-pe_%05d" \
        --cpus 4 \
        --ram 7; then

        # Check if error is due to no contigs
        if grep "has zero contigs" results/dragonflye.log; then
            touch results/contigs.fa
        else
            exit 1
        fi
    fi
    mv results/contigs.fa results/test-pe.fna
else
    # Shovill
    if ! shovill --R1 test-pe_R1.fastq.gz --R2 test-pe_R2.fastq.gz \
        --gsize 358000 \
        --outdir results \
        --assembler skesa --minlen 500 --mincov 2 --force --keepfiles --depth 0 --noreadcorr \
        --namefmt "test-pe_%05d" \
        --cpus 4 \
        --ram 7; then

        # Check if error is due to no contigs
        if grep "has zero contigs" results/shovill.log; then
            touch results/contigs.fa
        else
            exit 1
        fi
    fi
    mv results/contigs.fa results/test-pe.fna

    # Rename Graphs
    if [ -f "results/contigs.gfa" ]; then
        mv results/contigs.gfa results/skesa-unpolished.gfa
    elif [ -f "results/contigs.fastg" ]; then
        mv results/contigs.fastg results/skesa-unpolished.gfa
    elif [ -f "results/contigs.LastGraph" ]; then
        mv results/contigs.LastGraph results/skesa-unpolished.gfa
    fi

    if [ -f "results/flye-info.txt" ]; then
        mv results/flye-info.txt results/flye.log
    fi
fi

# Check quality of assembly
TOTAL_CONTIGS=`grep -c "^>" results/test-pe.fna || true`
if [ "${TOTAL_CONTIGS}" -gt 0 ]; then
    assembly-scan results/test-pe.fna --prefix test-pe > results/test-pe.tsv
    TOTAL_CONTIG_SIZE=$(cut -f 3 results/test-pe.tsv | tail -n 1)
    if [ "${TOTAL_CONTIG_SIZE}" -lt 100000 ]; then
        mv results/test-pe.fna results/test-pe-error.fna
        mv results/test-pe.tsv results/test-pe-error.tsv
        echo "test-pe assembled size (${TOTAL_CONTIG_SIZE} bp) is less than the minimum allowed genome
                size (100000 bp). If this is unexpected, please investigate test-pe to
                determine a cause (e.g. metagenomic, contaminants, etc...) for the poor assembly.
                Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                based analysis of test-pe will be discontinued." |             sed 's/^\s*//' > test-pe-assembly-error.txt
    fi
else
    mv results/test-pe.fna results/test-pe-error.fna
    echo "test-pe assembled successfully, but 0 contigs were formed. Please investigate
            test-pe to determine a cause (e.g. metagenomic, contaminants, etc...) for this
            outcome. Further assembly-based analysis of test-pe will be discontinued." |         sed 's/^\s*//' > test-pe-assembly-error.txt
fi

# Cleanup and compress
if [ "false" == "false" ]; then
    # Remove intermediate files
    rm -rfv results/shovill.bam* \
            results/shovill-se.bam* \
            results/flash.extendedFrags* \
            results/flash.notCombined* \
            results/skesa.fasta* \
            results/*.fq.gz \
            results/00*.gfa \
            results/pilon_polish* \
            results/flye/ \
            results/flye.fasta* \
            results/raven/ \
            results/raven.fasta* \
            results/raven.cereal \
            results/miniasm/ \
            results/miniasm.fasta* \
            results/spades/ \
            results/spades.fasta* \
            results/megahit/ \
            results/megahit.fasta* \
            results/velvet.fasta* \
            results/velvet/
fi

if [[ false == "false" ]]; then
    # Compress based on matched extensions
    find results/ -type f |             grep -E "\.fna$|\.fasta$|\.fa$|\.gfa$" |             xargs -I {} pigz -n --best -p 4 {}
fi
find results -maxdepth 1 -name "*.log" | xargs -I {} mv {} ./

# Capture versions
if [[ "$OSTYPE" == "darwin"* ]]; then

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ASSEMBLER:ASSEMBLER_MODULE":
    assembly-scan: $(echo $(assembly-scan --version 2>&1) | sed 's/assembly-scan //')
    bwa: $(echo $(bwa 2>&1) | sed 's/^.*Version: //;s/ .*$//')
    flash: $(echo $(flash --version 2>&1) | sed 's/^.*FLASH v//;s/ .*$//')
    megahit: $(echo $(megahit --version 2>&1) | sed 's/MEGAHIT v//')
    miniasm: $(echo $(miniasm -V))
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    pilon: $(echo $(pilon --version 2>&1) | sed 's/^.*Pilon version //;s/ .*$//')
    racon: $(echo $(racon --version 2>&1) | sed 's/v//')
    samclip: $(echo $(samclip --version 2>&1) | sed 's/samclip //')
    samtools: $(echo $(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*$//')
    shovill: $(echo $(shovill --version 2>&1) | sed 's/shovill //')
    shovill-se: $(echo $(shovill-se --version 2>&1) | sed 's/shovill-se //')
    skesa: $(echo $(skesa --version 2>&1) | sed 's/^.*SKESA //;s/ .*$//')
    spades.py: $(echo $(spades.py --version 2>&1) | sed 's/SPAdes genome assembler v//')
    velvetg: $(echo $(velvetg 2>&1) | sed 's/^.*Version //;s/ .*$//')
    velveth: $(echo $(velveth 2>&1) | sed 's/^.*Version //;s/ .*$//')
    unicycler: $(echo $(unicycler --version 2>&1) | sed 's/^.*Unicycler v//;s/ .*$//')
END_VERSIONS

else

cat <<-END_VERSIONS > versions.yml
"BACTOPIA:ASSEMBLER:ASSEMBLER_MODULE":
    any2fasta: $(echo $(any2fasta -v 2>&1) | sed 's/any2fasta //')
    assembly-scan: $(echo $(assembly-scan --version 2>&1) | sed 's/assembly-scan //')
    bwa: $(echo $(bwa 2>&1) | sed 's/^.*Version: //;s/ .*$//')
    flash: $(echo $(flash --version 2>&1) | sed 's/^.*FLASH v//;s/ .*$//')
    flye: $(echo $(flye --version))
    medaka: $(echo $(medaka --version 2>&1) | sed 's/medaka //')
    megahit: $(echo $(megahit --version 2>&1) | sed 's/MEGAHIT v//')
    miniasm: $(echo $(miniasm -V))
    minimap2: $(echo $(minimap2 --version))
    nanoq: $(echo $(nanoq --version 2>&1) | sed 's/nanoq //')
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    pilon: $(echo $(pilon --version 2>&1) | sed 's/^.*Pilon version //;s/ .*$//')
    racon: $(echo $(racon --version 2>&1) | sed 's/v//')
    rasusa: $(echo $(rasusa --version 2>&1) | sed 's/rasusa //')
    raven: $(echo $(raven --version))
    samclip: $(echo $(samclip --version 2>&1) | sed 's/samclip //')
    samtools: $(echo $(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*$//')
    shovill: $(echo $(shovill --version 2>&1) | sed 's/shovill //')
    shovill-se: $(echo $(shovill-se --version 2>&1) | sed 's/shovill-se //')
    skesa: $(echo $(skesa --version 2>&1) | sed 's/^.*SKESA //;s/ .*$//')
    spades.py: $(echo $(spades.py --version 2>&1) | sed 's/SPAdes genome assembler v//')
    velvetg: $(echo $(velvetg 2>&1) | sed 's/^.*Version //;s/ .*$//')
    velveth: $(echo $(velveth 2>&1) | sed 's/^.*Version //;s/ .*$//')
    unicycler: $(echo $(unicycler --version 2>&1) | sed 's/^.*Unicycler v//;s/ .*$//')
END_VERSIONS

fi

#!/bin/bash -ue
MERGED="multiple-read-sets-merged.txt"
mkdir -p fastqs
mkdir -p extra

if [ "hybrid" == "paired-end" ]; then
    # Paired-End Reads
    cp -L 001-r1 fastqs/test-hybrid_R1.fastq.gz
    cp -L 001-r2 fastqs/test-hybrid_R2.fastq.gz
    touch extra/empty.fna.gz
elif [ "hybrid" == "single-end" ]; then
    # Single-End Reads
    cp -L 001-r1 fastqs/test-hybrid.fastq.gz
    touch extra/empty.fna.gz
elif [ "hybrid" == "ont" ]; then
    # Nanopore reads
    cp -L 001-r1 fastqs/test-hybrid.fastq.gz
    touch extra/empty.fna.gz
elif  [ "hybrid" == "hybrid" ] || [ "hybrid" == "short_polish" ]; then 
    # Paired-End Reads
    cp -L 001-r1 fastqs/test-hybrid_R1.fastq.gz
    cp -L 001-r2 fastqs/test-hybrid_R2.fastq.gz
    cp -L ERR3772599.fastq.gz extra/test-hybrid.fastq.gz
elif [ "hybrid" == "merge-pe" ] || [ "hybrid" == "hybrid-merge-pe" ]; then 
    # Merge Paired-End Reads
    echo "This sample had reads merged." > ${MERGED}
    echo "R1:" >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"	"$9}' >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/test-hybrid_R1.fastq.gz
    echo "Merged R1:" >> ${MERGED}
    ls -l fastqs/test-hybrid_R1.fastq.gz | awk '{print $5"	"$9}' >> ${MERGED}

    echo "R2:" >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"	"$9}' >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/test-hybrid_R2.fastq.gz
    echo "Merged R2:" >> ${MERGED}
    ls -l fastqs/test-hybrid_R2.fastq.gz | awk '{print $5"	"$9}' >> ${MERGED}

    if [ "hybrid" == "hybrid-merge-pe" ]; then
        cp -L ERR3772599.fastq.gz extra/test-hybrid.fastq.gz
    else
        touch extra/empty.fna.gz
    fi
elif [ "hybrid" == "merge-se" ]; then 
    # Merge Single-End Reads
    echo "This sample had reads merged." > ${MERGED}
    echo "SE:" >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"	"$9}' >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/test-hybrid.fastq.gz
    echo "Merged SE:" >> ${MERGED}
    ls -l fastqs/test-hybrid.fastq.gz | awk '{print $5"	"$9}' >> ${MERGED}

    touch extra/empty.fna.gz
elif [ "hybrid" == "sra_accession" ] || [ "hybrid" == "sra_accession_ont" ]; then
    if [ "1" == "3" ]; then
        echo "Unable to download test-hybrid from both SRA and ENA 3 times. This may or may 
            not be a temporary connection issue. Rather than stop the whole Bactopia run, 
            further analysis of test-hybrid will be discontinued." | \
        sed 's/^\s*//' > test-hybrid-fastq-download-error.txt
        exit
    else
        # Download accession from ENA/SRA
        fastq-dl \
            --accession test-hybrid \
            --provider SRA \
            --cpus 4 \
            --outdir fastqs/ \
            --prefix test-hybrid \
            --group-by-experiment
        touch extra/empty.fna.gz
    fi 
elif [ "false" == "true" ]; then
    if [ "hybrid" == "assembly_accession" ]; then
        if [ "1" == "3" ]; then
            touch extra/empty.fna.gz
            echo "Unable to download test-hybrid from NCBI Assembly 3 times. This may or may
                not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                further analysis of test-hybrid will be discontinued." | \
            sed 's/^\s*//' > test-hybrid-assembly-download-error.txt
            exit
        else
            # Verify Assembly accession
            check-assembly-accession.py test-hybrid > accession.txt 2> check-assembly-accession.txt

            if [ -s "accession.txt" ]; then
                # Download from NCBI assembly and simulate reads
                mkdir fasta/
                ncbi-genome-download bacteria -o ./ -F fasta -p 4 \
                                            -s null -A accession.txt -r 50 
                find . -name "*test-hybrid*.fna.gz" | xargs -I {} mv {} fasta/
                rename 's/(GC[AF]_\d+).*/$1.fna.gz/' fasta/*
                gzip -cd fasta/test-hybrid.fna.gz > test-hybrid-art.fna
                rm check-assembly-accession.txt
            else
                mv check-assembly-accession.txt test-hybrid-assembly-accession-error.txt
                exit
            fi
        fi
    elif [ "hybrid" == "assembly" ]; then
        if [ "true" == "true" ]; then
            gzip -cd ERR3772599.fastq.gz > test-hybrid-art.fna
        else 
            cat ERR3772599.fastq.gz > test-hybrid-art.fna
        fi
    fi

    # Simulate reads from assembly, reads are 250bp without errors
    art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov 150 -ir 0 -ir2 0 -dr 0 -dr2 0 -rs 42                        -na -qL 33 -qU 40 -o test-hybrid_R --id test-hybrid -i test-hybrid-art.fna

    mv test-hybrid_R1.fq fastqs/test-hybrid_R1.fastq
    mv test-hybrid_R2.fq fastqs/test-hybrid_R2.fastq
    pigz -p 4 --fast fastqs/*.fastq
    cp test-hybrid-art.fna extra/test-hybrid.fna
    pigz -p 4 --best extra/test-hybrid.fna
fi

# Validate input FASTQs
if [ "false" == "false" ]; then
    ERROR=0
    # Check paired-end reads have same read counts
    OPTS="--sample test-hybrid --min_basepairs 2241820 --min_reads 7472 --min_proportion 0.5 --runtype hybrid"
    if [ -f  "fastqs/test-hybrid_R2.fastq.gz" ]; then
        # Paired-end
        gzip -cd fastqs/test-hybrid_R1.fastq.gz | fastq-scan > r1.json
        gzip -cd fastqs/test-hybrid_R2.fastq.gz | fastq-scan > r2.json
        if ! reformat.sh in1=fastqs/test-hybrid_R1.fastq.gz in2=fastqs/test-hybrid_R2.fastq.gz qin=auto out=/dev/null 2> test-hybrid-paired-end-error.txt; then
            ERROR=1
            echo "test-hybrid FASTQs contains an error. Please check the input FASTQs.
                Further analysis is discontinued." | \
            sed 's/^\s*//' >> test-hybrid-paired-end-error.txt
        else
            rm -f test-hybrid-paired-end-error.txt
        fi

        if ! check-fastqs.py --fq1 r1.json --fq2 r2.json ${OPTS}; then
            ERROR=1
        fi
        rm r1.json r2.json
    else
        # Single-end
        gzip -cd fastqs/test-hybrid.fastq.gz | fastq-scan > r1.json
        if ! check-fastqs.py --fq1 r1.json ${OPTS}; then
            ERROR=1
        fi
        rm r1.json
    fi

    # Failed validations so, let's keep them from continuing
    if [ "${ERROR}" -eq "1" ]; then
        mv fastqs/ failed-tests-fastqs/
    fi
fi

# Dump meta values to a TSV
echo "sample<TAB>runtype<TAB>original_runtype<TAB>species<TAB>genome_size" | sed 's/<TAB>/	/g' > test-hybrid-meta.tsv
echo "test-hybrid<TAB>hybrid<TAB>hybrid<TAB>null<TAB>358000" | sed 's/<TAB>/	/g' >> test-hybrid-meta.tsv

# Capture versions
cat <<-END_VERSIONS > versions.yml
"BACTOPIA:GATHER:GATHER_MODULE":
    art: $(echo $(art_illumina --help 2>&1) | sed 's/^.*Version //;s/ .*$//')
    fastq-dl: $(echo $(fastq-dl --version 2>&1) | sed 's/fastq-dl, version //')
    fastq-scan: $(echo $(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
    mash: $(echo $(mash 2>&1) | sed 's/^.*Mash version //;s/ .*$//')
    ncbi-genome-download: $(echo $(ncbi-genome-download --version 2>&1))
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
END_VERSIONS

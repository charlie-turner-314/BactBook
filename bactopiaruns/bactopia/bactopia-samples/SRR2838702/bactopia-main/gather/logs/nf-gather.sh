#!/bin/bash -ue
MERGED="multiple-read-sets-merged.txt"
mkdir -p fastqs
mkdir -p extra

if [ "paired-end" == "paired-end" ]; then
    # Paired-End Reads
    cp -L 001-r1 fastqs/SRR2838702_R1.fastq.gz
    cp -L 001-r2 fastqs/SRR2838702_R2.fastq.gz
    touch extra/empty.fna.gz
elif [ "paired-end" == "single-end" ]; then
    # Single-End Reads
    cp -L 001-r1 fastqs/SRR2838702.fastq.gz
    touch extra/empty.fna.gz
elif [ "paired-end" == "ont" ]; then
    # Nanopore reads
    cp -L 001-r1 fastqs/SRR2838702.fastq.gz
    touch extra/empty.fna.gz
elif  [ "paired-end" == "hybrid" ] || [ "paired-end" == "short_polish" ]; then 
    # Paired-End Reads
    cp -L 001-r1 fastqs/SRR2838702_R1.fastq.gz
    cp -L 001-r2 fastqs/SRR2838702_R2.fastq.gz
    cp -L EMPTY_EXTRA extra/SRR2838702.fastq.gz
elif [ "paired-end" == "merge-pe" ] || [ "paired-end" == "hybrid-merge-pe" ]; then 
    # Merge Paired-End Reads
    echo "This sample had reads merged." > ${MERGED}
    echo "R1:" >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"	"$9}' >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/SRR2838702_R1.fastq.gz
    echo "Merged R1:" >> ${MERGED}
    ls -l fastqs/SRR2838702_R1.fastq.gz | awk '{print $5"	"$9}' >> ${MERGED}

    echo "R2:" >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"	"$9}' >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/SRR2838702_R2.fastq.gz
    echo "Merged R2:" >> ${MERGED}
    ls -l fastqs/SRR2838702_R2.fastq.gz | awk '{print $5"	"$9}' >> ${MERGED}

    if [ "paired-end" == "hybrid-merge-pe" ]; then
        cp -L EMPTY_EXTRA extra/SRR2838702.fastq.gz
    else
        touch extra/empty.fna.gz
    fi
elif [ "paired-end" == "merge-se" ]; then 
    # Merge Single-End Reads
    echo "This sample had reads merged." > ${MERGED}
    echo "SE:" >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"	"$9}' >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/SRR2838702.fastq.gz
    echo "Merged SE:" >> ${MERGED}
    ls -l fastqs/SRR2838702.fastq.gz | awk '{print $5"	"$9}' >> ${MERGED}

    touch extra/empty.fna.gz
elif [ "paired-end" == "sra_accession" ] || [ "paired-end" == "sra_accession_ont" ]; then
    if [ "1" == "3" ]; then
        echo "Unable to download SRR2838702 from both SRA and ENA 3 times. This may or may 
            not be a temporary connection issue. Rather than stop the whole Bactopia run, 
            further analysis of SRR2838702 will be discontinued." | \
        sed 's/^\s*//' > SRR2838702-fastq-download-error.txt
        exit
    else
        # Download accession from ENA/SRA
        fastq-dl \
            --accession SRR2838702 \
            --provider SRA \
            --cpus 2 \
            --outdir fastqs/ \
            --prefix SRR2838702 \
            --group-by-experiment
        touch extra/empty.fna.gz
    fi 
elif [ "false" == "true" ]; then
    if [ "paired-end" == "assembly_accession" ]; then
        if [ "1" == "3" ]; then
            touch extra/empty.fna.gz
            echo "Unable to download SRR2838702 from NCBI Assembly 3 times. This may or may
                not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                further analysis of SRR2838702 will be discontinued." | \
            sed 's/^\s*//' > SRR2838702-assembly-download-error.txt
            exit
        else
            # Verify Assembly accession
            check-assembly-accession.py SRR2838702 > accession.txt 2> check-assembly-accession.txt

            if [ -s "accession.txt" ]; then
                # Download from NCBI assembly and simulate reads
                mkdir fasta/
                ncbi-genome-download bacteria -o ./ -F fasta -p 2 \
                                            -s null -A accession.txt -r 50 
                find . -name "*SRR2838702*.fna.gz" | xargs -I {} mv {} fasta/
                rename 's/(GC[AF]_\d+).*/$1.fna.gz/' fasta/*
                gzip -cd fasta/SRR2838702.fna.gz > SRR2838702-art.fna
                rm check-assembly-accession.txt
            else
                mv check-assembly-accession.txt SRR2838702-assembly-accession-error.txt
                exit
            fi
        fi
    elif [ "paired-end" == "assembly" ]; then
        if [ "false" == "true" ]; then
            gzip -cd EMPTY_EXTRA > SRR2838702-art.fna
        else 
            cat EMPTY_EXTRA > SRR2838702-art.fna
        fi
    fi

    # Simulate reads from assembly, reads are 250bp without errors
    art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov 150 -ir 0 -ir2 0 -dr 0 -dr2 0 -rs 42                        -na -qL 33 -qU 40 -o SRR2838702_R --id SRR2838702 -i SRR2838702-art.fna

    mv SRR2838702_R1.fq fastqs/SRR2838702_R1.fastq
    mv SRR2838702_R2.fq fastqs/SRR2838702_R2.fastq
    pigz -p 2 --fast fastqs/*.fastq
    cp SRR2838702-art.fna extra/SRR2838702.fna
    pigz -p 2 --best extra/SRR2838702.fna
fi

# Validate input FASTQs
if [ "false" == "false" ]; then
    ERROR=0
    # Check paired-end reads have same read counts
    OPTS="--sample SRR2838702 --min_basepairs 2241820 --min_reads 7472 --min_proportion 0.5 --runtype paired-end"
    if [ -f  "fastqs/SRR2838702_R2.fastq.gz" ]; then
        # Paired-end
        gzip -cd fastqs/SRR2838702_R1.fastq.gz | fastq-scan > r1.json
        gzip -cd fastqs/SRR2838702_R2.fastq.gz | fastq-scan > r2.json
        if ! reformat.sh in1=fastqs/SRR2838702_R1.fastq.gz in2=fastqs/SRR2838702_R2.fastq.gz qin=auto out=/dev/null 2> SRR2838702-paired-end-error.txt; then
            ERROR=1
            echo "SRR2838702 FASTQs contains an error. Please check the input FASTQs.
                Further analysis is discontinued." | \
            sed 's/^\s*//' >> SRR2838702-paired-end-error.txt
        else
            rm -f SRR2838702-paired-end-error.txt
        fi

        if ! check-fastqs.py --fq1 r1.json --fq2 r2.json ${OPTS}; then
            ERROR=1
        fi
        rm r1.json r2.json
    else
        # Single-end
        gzip -cd fastqs/SRR2838702.fastq.gz | fastq-scan > r1.json
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
echo "sample<TAB>runtype<TAB>original_runtype<TAB>species<TAB>genome_size" | sed 's/<TAB>/	/g' > SRR2838702-meta.tsv
echo "SRR2838702<TAB>paired-end<TAB>paired-end<TAB>null<TAB>358242" | sed 's/<TAB>/	/g' >> SRR2838702-meta.tsv

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

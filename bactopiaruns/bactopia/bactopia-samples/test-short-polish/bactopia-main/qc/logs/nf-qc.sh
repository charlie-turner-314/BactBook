#!/bin/bash -ue
mkdir -p results
touch results/.short_polish
ERROR=0
MIN_COVERAGE=$(( 10*358000 ))
TOTAL_BP=$(( 100*358000 ))

if [[ "false" == "true" ]]; then
    echo "Sequence QC was skipped for test-short-polish" > results/test-short-polish-qc-skipped.txt
    if [ "false" == "false" ]; then
        # Paired-End Reads
        cp test-short-polish_R1.fastq.gz results/test-short-polish_R1.fastq.gz
        cp test-short-polish_R2.fastq.gz results/test-short-polish_R2.fastq.gz
    else
        # Single-End Reads
        cp test-short-polish_R1.fastq.gz results/test-short-polish.fastq.gz
    fi
else
    if [[ "short_polish" == "ont" ]]; then
        if [[ "false" == "true" ]]; then
            # Remove Adapters
            porechop --input test-short-polish_R1.fastq.gz                      --format fastq                     --threads 4 > adapter-r1.fq

            # Quality filter
            nanoq --min-len 1000                     --min-qual 0                     --input adapter-r1.fq 1> filt-r1.fq
        else 
            # Quality filter
            nanoq --min-len 1000                     --min-qual 0                     --input  test-short-polish_R1.fastq.gz 1> filt-r1.fq
        fi
    elif [[ "false" == "true" ]]; then
        # Use BBMap for cleaning reads
        # Illumina Reads
        # Validate paired-end reads if necessary
        if [[ "false" == "false" ]]; then
            # Make sure paired-end reads have matching IDs
            repair.sh                     in=test-short-polish_R1.fastq.gz                     in2=test-short-polish_R2.fastq.gz                     out=repair-r1.fq                     out2=repair-r2.fq                     outs=repair-singles.fq                     ain=f

            if [ ! -s repair-r1.fq ]; then
                ERROR=1
                echo "After validating read pairs, test-short-polish FASTQs are empty. Please check the input FASTQs.
                    Further analysis is discontinued." |                     sed 's/^\s*//' >> test-short-polish-paired-match-error.txt
            fi
        else
            gunzip -c test-short-polish_R1.fastq.gz > repair-r1.fq 
        fi

        if [ "${ERROR}" -eq "0" ]; then
            # Remove Adapters
            bbduk.sh -Xmx8160437862                     in=repair-r1.fq out=adapter-r1.fq in2=repair-r2.fq out2=adapter-r2.fq                     ref=adapters                     k=23                     ktrim=r                     mink=11                     hdist=1                     tpe=t                     tbo=t                     threads=4                     ftm=5                     qin=auto ordered=t 

            if [ ! -s adapter-r1.fq ]; then
                ERROR=1
                echo "After adapter removal, test-short-polish FASTQs are empty. Please check the input FASTQs.
                    Further analysis is discontinued." |                     sed 's/^\s*//' >> test-short-polish-adapter-qc-error.txt
            fi
        fi

        if [ "${ERROR}" -eq "0" ]; then
            # Remove PhiX
            bbduk.sh -Xmx8160437862                     in=adapter-r1.fq out=phix-r1.fq in2=adapter-r2.fq out2=phix-r2.fq                     ref=phix                     k=31                     hdist=1                     tpe=t                     tbo=t                     qtrim=rl                     trimq=6                     minlength=35                     minavgquality=10                     qin=auto qout=33                     tossjunk=t                     threads=4                     ordered=t 

            if [ ! -s phix-r1.fq ]; then
                ERROR=1
                echo "After PhiX removal, test-short-polish FASTQs are empty. Please check the input FASTQs.
                    Further analysis is discontinued." |                     sed 's/^\s*//' >> test-short-polish-phix-qc-error.txt
            fi
        fi
    else
        # QC with fastp
        mkdir -p results/summary/
        fastp                 --in1 test-short-polish_R1.fastq.gz --out1 filt-r1.fq --in2 test-short-polish_R2.fastq.gz --out2 filt-r2.fq --detect_adapter_for_pe                 --thread 4                 --json results/summary/test-short-polish.fastp.json                 --html results/summary/test-short-polish.fastp.html  2> test-short-polish-fastp.log
    fi

    # Error Correction
    if [ "${ERROR}" -eq "0" ]; then
        if [[ "false" == "true" ]]; then
            if [ "false" == "false" ] && [ "358000" -gt "0" ]; then
                lighter -od . -r phix-r1.fq -r phix-r2.fq -K 31 358000 -maxcor 1 -zlib 0 -t 4
                mv phix-r1.cor.fq filt-r1.fq
                if [[ "false" == "false" ]]; then
                    mv phix-r2.cor.fq filt-r2.fq
                fi
            else
                echo "Skipping error correction"
                ln -s phix-r1.fq filt-r1.fq
                if [ "false" == "false" ]; then
                    ln -s phix-r2.fq filt-r2.fq
                fi
            fi
        elif [[ "short_polish" == "ont" ]]; then
            echo "Skipping error correction. Have a recommended ONT error corrector? Let me know!"
        else
            echo "Skipping error correction since fastp was used"
        fi
    fi

    # Reduce Coverage
    if [ "${ERROR}" -eq "0" ]; then
        if (( ${TOTAL_BP} > 0 )); then
            if [[ "short_polish" == "ont" ]]; then
                rasusa -i filt-r1.fq                         -c 100                         -g 358000                         -s 42 1> subsample-r1.fq
            else
                if [ -f filt-r1.fq ]; then
                    reformat.sh -Xmx8160437862                             in=filt-r1.fq out=subsample-r1.fq in2=filt-r2.fq out2=subsample-r2.fq                             samplebasestarget=${TOTAL_BP}                             sampleseed=42                             overwrite=t
                fi
            fi
        else
            echo "Skipping coverage reduction"
            ln -s filt-r1.fq subsample-r1.fq
            if [ "false" == "false" ]; then
                ln -s filt-r2.fq subsample-r2.fq
            fi
        fi
    fi

    # Compress
    if [ "${ERROR}" -eq "0" ]; then
        if [ "false" == "false" ]; then
            pigz -p 4 -c -n subsample-r1.fq > results/test-short-polish_R1.fastq.gz
            pigz -p 4 -c -n subsample-r2.fq > results/test-short-polish_R2.fastq.gz
        else
            pigz -p 4 -c -n subsample-r1.fq > results/test-short-polish.fastq.gz
        fi

        if [ "false" == "false" ]; then
            # Remove remaining intermediate FASTQ files
            rm *.fq
        fi
    fi
fi

# Quality stats before and after QC
if [ "${ERROR}" -eq "0" ]; then
    mkdir -p results/summary/
    # fastq-scan
    if [[ "false" == "false" ]]; then
        # Paired-End Reads
        gzip -cd test-short-polish_R1.fastq.gz | fastq-scan -g 358000 > results/summary/test-short-polish_R1-original.json
        gzip -cd test-short-polish_R2.fastq.gz | fastq-scan -g 358000 > results/summary/test-short-polish_R2-original.json
        gzip -cd results/test-short-polish_R1.fastq.gz | fastq-scan -g 358000 > results/summary/test-short-polish_R1-final.json
        gzip -cd results/test-short-polish_R2.fastq.gz | fastq-scan -g 358000 > results/summary/test-short-polish_R2-final.json
    else
        # Single-End Reads
        gzip -cd test-short-polish_R1.fastq.gz | fastq-scan -g 358000 > results/summary/test-short-polish-original.json
        gzip -cd results/test-short-polish.fastq.gz | fastq-scan -g 358000 > results/summary/test-short-polish-final.json
    fi

    # FastQC and NanoPlot
    if [[ "false" == "false" ]]; then
        if [[ "short_polish" == "ont" ]]; then
            mkdir results/summary/test-short-polish-original results/summary/test-short-polish-final
            NanoPlot                      --threads 4                     --fastq test-short-polish_R1.fastq.gz                     --outdir results/summary/test-short-polish-original/                     --prefix test-short-polish-original_
            cp results/summary/test-short-polish-original/test-short-polish-original_NanoPlot-report.html results/summary/test-short-polish-original_NanoPlot-report.html
            tar -cvf - results/summary/test-short-polish-original/ | pigz --best -p 4 > results/summary/test-short-polish-original_NanoPlot.tar.gz

            NanoPlot                      --threads 4                     --fastq results/test-short-polish.fastq.gz                     --outdir results/summary/test-short-polish-final/                     --prefix test-short-polish-final_
            cp results/summary/test-short-polish-final/test-short-polish-final_NanoPlot-report.html results/summary/test-short-polish-final_NanoPlot-report.html
            tar -cvf - results/summary/test-short-polish-final/ | pigz --best -p 4 > results/summary/test-short-polish-final_NanoPlot.tar.gz
            rm -rf results/summary/test-short-polish-original/ results/summary/test-short-polish-final/
        else
            if [ "false" == "false" ]; then
                # Paired-End Reads
                ln -s test-short-polish_R1.fastq.gz test-short-polish_R1-original.fastq.gz
                ln -s test-short-polish_R2.fastq.gz test-short-polish_R2-original.fastq.gz
                ln -s results/test-short-polish_R1.fastq.gz test-short-polish_R1-final.fastq.gz
                ln -s results/test-short-polish_R2.fastq.gz test-short-polish_R2-final.fastq.gz
                fastqc --noextract -f fastq -t 4 test-short-polish_R1-original.fastq.gz test-short-polish_R2-original.fastq.gz test-short-polish_R1-final.fastq.gz test-short-polish_R2-final.fastq.gz
            else
                # Single-End Reads
                ln -s test-short-polish_R1.fastq.gz test-short-polish-original.fastq.gz
                ln -s results/test-short-polish.fastq.gz test-short-polish-final.fastq.gz
                fastqc --noextract -f fastq -t 4 test-short-polish-original.fastq.gz test-short-polish-final.fastq.gz
            fi
            mv *_fastqc.html *_fastqc.zip results/summary/
        fi
    fi
fi

# Final QC check
if [ "false" == "false" ]; then
    # Only check for errors if we haven't already found them
    if [ "${ERROR}" -eq "0" ]; then
        gzip -cd results/*.fastq.gz | fastq-scan -g 358000 > temp.json
        FINAL_BP=$(grep "total_bp" temp.json | sed -r 's/.*:[ ]*([0-9]+),/\1/')
        rm temp.json

        if [ ${FINAL_BP} -lt ${MIN_COVERAGE} ]; then
            ERROR=1
            echo "After QC, test-short-polish FASTQ(s) contain ${FINAL_BP} total basepairs. This does
                    not exceed the required minimum ${MIN_COVERAGE} bp (10x coverage). Further analysis 
                    is discontinued." |                 sed 's/^\s*//' > test-short-polish-low-sequence-coverage-error.txt
            ERROR=2
        fi

        # Check paired-end reads have same read counts
        OPTS="--sample test-short-polish --min_basepairs 2241820 --min_reads 7472 --min_proportion 0.5 --runtype short_polish"
        if [ -f  "results/test-short-polish_R2.fastq.gz" ]; then
            # Paired-end
            gzip -cd results/test-short-polish_R1.fastq.gz | fastq-scan > r1.json
            gzip -cd results/test-short-polish_R2.fastq.gz | fastq-scan > r2.json
            if ! reformat.sh in1=results/test-short-polish_R1.fastq.gz in2=results/test-short-polish_R2.fastq.gz qin=auto out=/dev/null 2> test-short-polish-paired-end-error.txt; then
                ERROR=2
                echo "test-short-polish FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." |                     sed 's/^\s*//' >> test-short-polish-paired-end-error.txt
            else
                rm -f test-short-polish-paired-end-error.txt
            fi
            if ! check-fastqs.py --fq1 r1.json --fq2 r2.json ${OPTS}; then
                ERROR=2
            fi
            rm r1.json r2.json
        else
            # Single-end
            gzip -cd results/test-short-polish.fastq.gz | fastq-scan > r1.json
            if ! check-fastqs.py --fq1 r1.json ${OPTS}; then
                ERROR=2
            fi
            rm r1.json
        fi
    fi
fi

if [ "false" == "true" ]; then
    touch results/reads-simulated-from-assembly.txt
fi

if [ "${ERROR}" -eq "1" ]; then
    if [ "false" == "false" ]; then
        cp test-short-polish_R1.fastq.gz results/test-short-polish_R1.error-fastq.gz
        cp test-short-polish_R2.fastq.gz results/test-short-polish_R2.error-fastq.gz
        if [ ! -s repair-singles.fq ]; then
            pigz -p 4 -c -n repair-singles.fq > results/test-short-polish.error-fastq.gz
        fi
    else
        cp test-short-polish_R1.fastq.gz results/test-short-polish.error-fastq.gz
    fi
elif [ "${ERROR}" -eq "2" ]; then
    if [ "false" == "false" ]; then
        if [ -f results/test-short-polish_R1.fastq.gz ]; then
            mv results/test-short-polish_R1.fastq.gz results/test-short-polish_R1.error-fastq.gz
            mv results/test-short-polish_R2.fastq.gz results/test-short-polish_R2.error-fastq.gz

            if [ -s repair-singles.fq ]; then
                pigz -p 4 -c -n repair-singles.fq > results/test-short-polish.error-fastq.gz
            fi
        fi
    else
        if [ -f results/test-short-polish.fastq.gz ]; then
            mv results/test-short-polish.fastq.gz results/test-short-polish.error-fastq.gz
        fi
    fi
fi

# Capture versions
if [[ "false" == "false" ]]; then
cat <<-END_VERSIONS > versions.yml
"BACTOPIA:QC:QC_MODULE":
    bbduk: $(echo $(bbduk.sh --version 2>&1) | sed 's/^.*BBMap version //;s/ .*$//')
    fastp: $(echo $(fastp --version 2>&1) | sed -e "s/fastp //g")
    fastqc: $(echo $(fastqc --version 2>&1) | sed 's/^.*FastQC v//')
    fastq-scan: $(echo $(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
    lighter: $(echo $(lighter -v 2>&1) | sed 's/Lighter v//')
    nanoplot: $(echo $(NanoPlot -v 2>&1) | sed 's/.*NanoPlot //')
    nanoq: $(echo $(nanoq --version 2>&1) | sed 's/nanoq //')
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    porechop: $(echo $(porechop --version 2>&1))
    rasusa: $(echo $(rasusa --version 2>&1) | sed 's/rasusa //')
END_VERSIONS
else
cat <<-END_VERSIONS > versions.yml
"BACTOPIA:QC:QC_MODULE":
    bbduk: $(echo $(bbduk.sh --version 2>&1) | sed 's/^.*BBMap version //;s/ .*$//')
    fastp: $(echo $(fastp --version 2>&1) | sed -e "s/fastp //g")
    fastq-scan: $(echo $(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
    lighter: $(echo $(lighter -v 2>&1) | sed 's/Lighter v//')
    nanoq: $(echo $(nanoq --version 2>&1) | sed 's/nanoq //')
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    porechop: $(echo $(porechop --version 2>&1))
    rasusa: $(echo $(rasusa --version 2>&1) | sed 's/rasusa //')
END_VERSIONS
fi

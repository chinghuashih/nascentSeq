#!/bin/bash
#SBATCH -p standard  -o proseq.mapping.bowtie2.log -t 8:00:00
#SBATCH -c 8 --mem=128G

module load fastqc
module load fastp
module load bowtie2
module load samtools
#module load deeptools/2.5.3
module load deeptools
module load umi-tools/b1

###########################################################
# This is a pipeline script for handling paired end       #
# PRO-seq data with UMIs on both ends of the read.        #
# Run this script in a directory that has one folder      #
# named "fastq" which contains the data.                  #
# Fastq files must have identical names other than        #
# ending in _R1.fastq and _R2.fastq.                      #
###########################################################

## Parameters
THREADS=8 # Threads to use for multithreaded applications
UMI_LEN=6 # Length of UMI in basepairs
SPIKEIN="Y"

## UMI Flags (set to Y or N as appropriate)
FIVEP_UMI="Y"  # Is there a UMI on the 5' end of the read?
THREEP_UMI="Y" # Is there a UMI on the 3' end of the read?

## Adaptor sequences to clip. Default = Tru-Seq small RNA
#ADAPTOR_1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" 
#ADAPTOR_2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

ADAPTOR_1="TGGAATTCTCGGGTGCCAAGG"
ADAPTOR_2="GATCGTCGGACTGTAGAACTCTGAAC"

## Genomes. Fill in paths.
GENOME_EXP="/scratch/cshih8/references/hg38/bowtie2_index/default/hg38"
RDNA="/scratch/cshih8/references/human_rDNA/bowtie2_index/default/human_rDNA"

GENOME_SPIKE="/scratch/cshih8/references/drosophila/dm6/dm6" ## USE REPEAT MASKED VERSION!!
SPIKE_PREFIX="dm6" ## This is the prefix you've used on your spike in chromosomes, ie >hg38chr1
#RDNA="/scratch/cshih8/references/hg38dm6/hg38dm6_rDNA"

## Mapq value for filtering multimappers
MAPQ=10
binSize=10

###############################################################
#                           PIPELINE                          #
###############################################################

mkdir -p logs/align
mkdir -p logs/deDup
mkdir -p logs/fastp
mkdir -p logs/fastqc
mkdir -p logs/rRNA
mkdir -p logs/spikeAlign
mkdir -p logs/spikedeDup

mkdir -p BAM
mkdir -p BAMdeDuped
mkdir -p spikeBAM
mkdir -p spikeBAMdeDuped

mkdir -p bigWig
mkdir -p info
mkdir -p trimmedFastq

# Unzipping if needed
echo "unzipping..."
for FILE in fastq/*
do
	if [[ "${FILE}" == *.gz ]]; then
		gunzip ${FILE} &
	fi
done
wait

# Removing extra info from filenames.
# This is general and works with files from Cornell BRC.
# If filenames are formatted differently, does nothing.
echo "renaming if needed..."
for FILE in $(ls fastq/)
do
	NEW=fastq/"$(echo "$FILE" | 
		sed 's/^[0-9]\+_[0-9]\+_[0-9]\+_[0-9A-Z]\+_//' | 
		sed 's/_[ATCG]\{6,8\}_/_/')"
	if [ ! -s "${NEW}" ]; then
		mv fastq/${FILE} ${NEW}
	fi
done

### Running fastqc on files
echo "running fastqc if needed..."
for FILE in fastq/*.fastq
do
	if [ ! -s logs/fastqc/"$(basename ${FILE/.fastq/_fastqc.zip})" ]; then
		fastqc ${FILE} -o logs/fastqc --quiet  &
	fi
done
wait

### Autodetecting paired end files
echo "detecting the type of libraries"
NUM=$(ls fastq | wc -l) 
NUM_REDUCED=$(ls fastq | sed 's/_R.*//' | uniq | wc -l)
if [[ ${NUM} == ${NUM_REDUCED} ]]; then
	PAIRED="N"
	echo "detected ""${NUM}"" single end fastq files."
else
	PAIRED="Y"
	echo "detected ""${NUM_REDUCED}"" paired end fastq files"
fi


if [[ ${PAIRED} == "N" ]]; then
	### Trimming adapters and filtering rRNA reads
	echo "trimming adapters and filtering rDNA reads..."

	## 3' UMI
	if [[ ${THREEP_UMI} == "Y" ]]; then
		for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
		do
			if [ ! -s trimmedFastq/${PAIR}.fastq ]; then
				(fastp \
					-i fastq/${PAIR}_R1.fastq \
					--adapter_sequence ${ADAPTOR_1} \
					--umi \
					--stdout \
					--umi_loc=read1 \
					--umi_len=${UMI_LEN} \
					--html logs/fastp/${PAIR}_fastp.html \
					-w ${THREADS} \
					--overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
				(bowtie2 \
					--fast-local \
					--un trimmedFastq/${PAIR}_R1.fastq \
					-U - \
					-x ${RDNA} \
					--threads ${THREADS} 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
				mv fastp.json logs/fastp/${PAIR}_fastp.json
			fi
		done

	# Branch for no UMI
	else
		for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
		do
			if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]; then
				(fastp \
					-i fastq/${PAIR}_R1.fastq \
					--adapter_sequence ${ADAPTOR_1} \
					--stdout \
					--html logs/fastp/${PAIR}_fastp.html \
					-w ${THREADS} \
					--overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
				(bowtie2 \
					--fast-local \
					--un trimmedFastq/${PAIR}_R1.fastq \
					--interleaved - \
					-x ${RDNA} \
					--threads ${THREADS} 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
				mv fastp.json logs/fastp/${PAIR}_fastp.json
			fi
		done
	fi

	### Cleaning up filenames in trimmedFastq
#	for FILE in trimmedFastq/*.fastq
#	do
#		if [ ! -s ${FILE/.fastq/_R1.fastq} ]; then
#			mv "${FILE}" ${FILE/.fastq/_R1.fastq}
#		fi
#	done

	### Aligning to spike in genome to get normalization factors
	for PAIR in $(ls trimmedFastq | sed 's/_R[1-2].*//' | uniq )
	do
		if [ ! -s "spikeBAM/${PAIR}.BAM" ]; then
			echo "aligning ${PAIR} to spike in genome"
			(bowtie2 \
				--local \
				--very-sensitive-local \
				--threads ${THREADS} \
				--no-unal \
				--no-mixed \
				--no-discordant \
				-x "${GENOME_SPIKE}" \
				-U trimmedFastq/${PAIR}_R1.fastq \
				2> logs/spikeAlign/MCF10A_CTRL_A_spikeAlign.log) |
			samtools view -hS -q ${MAPQ} |
			samtools view -b | 
			samtools sort -@ ${THREADS} -o spikeBAM/${PAIR}.BAM
			samtools index spikeBAM/${PAIR}.BAM

			#perl -n -e 'print $_ if (/^\@/ || /'${SPIKE_PREFIX}'/ ) ' |
		fi
	done

	### Aligning to experimental genome

	for PAIR in $(ls trimmedFastq | sed 's/_R[1-2].*//' | uniq )
	do
		if [ ! -s "BAM/${PAIR}.BAM" ]; then
			echo "aligning ${PAIR} to experimental genome"
			(bowtie2 \
				--local \
				--sensitive-local \
				--threads ${THREADS} \
				-x "${GENOME_EXP}" \
				-U "trimmedFastq/${PAIR}_R1.fastq" \
				2> logs/align/${PAIR}_align.log) |
			samtools view -bS -q ${MAPQ} |
			samtools sort -@ ${THREADS} -o BAM/${PAIR}.BAM
			samtools index BAM/${PAIR}.BAM
			echo ""
		fi
	done

	### deduplicating with UMIs
	module unload python
	module load python3

	for FILE in BAM/*.BAM
	do
		echo "dedup exp"
		
		if [ ! -s "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" ]; then
			(umi_tools dedup \
				-I "${FILE}" \
				--umi-separator=":" \
				-S "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
				)> "logs/deDup/$(basename ${FILE%.BAM}_deDup.log)" &&
			samtools index "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" 
		fi
	done

	# deduplicating with UMIs

	for FILE in spikeBAM/*.BAM
	do
		echo "dedup spike-in"
		if [ ! -s "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" ]; then
			(umi_tools dedup \
				  -I "${FILE}" \
				  --umi-separator=":" \
				  -S "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
				  )> "logs/spikedeDup/$(basename ${FILE%.BAM}_deDup.log)" &&
				  samtools index "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" 
		fi
	done

	### Generating infoTable

	echo "stat"
	if [ ! -s info/infoTable.tsv ]; then
		touch info/infoTable.tsv
		echo -e \
			Name'\t'\
			RawReads'\t'\
			NonDimerReads'\t'\
			%dimer'\t'\
			rRNAreads'\t'\
			%rRNA'\t'\
			passedFilters'\t'\
			bowtieConcordant'\t'\
			bowtieMulti'\t'\
			bowtieUnal'\t'\
			bowtieOverallMap%'\t'\
			bowtieConcordant%'\t'\
			bowtieMulti%'\t'\
			bowtieUnal%'\t'\
			uniqueMapped'\t'\
			uniqueMappedNondup'\t'\
			%PCRdups'\t'\
			uniqueMappedSpikein'\t'\
			uniqueMappedSpikeinNondup'\t'\
			spikeInPCRdups% >> info/infoTable.tsv

		for SAMPLE in $(ls BAM/*.BAM | sed 's/.BAM//' | sed 's/BAM\///' )
		do
			NAME=${SAMPLE}
			RAW_READS=$(cat logs/fastp/${SAMPLE}_fastp.log              |
					grep "total reads:"                         | head -n 1 | awk '{print $3}')
			TRIMMED_READS=$(cat logs/fastp/${SAMPLE}_fastp.log          |
					grep "total reads:"                         | tail -n 1 | awk '{print $3}')
			PER_DIMER=$(echo "(1-"${TRIMMED_READS}"/"${RAW_READS}")*100" | bc -l)%
			PASSED_FILTERS=$(cat logs/align/${SAMPLE}_align.log         |
					grep "reads; of these:$"                    | awk '{print $1}')
			RRNA=$(echo ${TRIMMED_READS}"-"${PASSED_FILTERS} | bc )
			PER_RRNA=$(echo ${RRNA}"/"${RAW_READS}"*100" | bc -l)%
			B_CONC=$(cat logs/align/${SAMPLE}_align.log                 |
					grep "aligned exactly 1 time$" | awk '{print $1}')
			B_MULTI=$(cat logs/align/${SAMPLE}_align.log                |
					grep "aligned >1 times$"       | awk '{print $1}')
			B_UNAL=$(cat logs/align/${SAMPLE}_align.log                 |
					grep "aligned 0 times$"        | awk '{print $1}')
			B_OAP=$(cat logs/align/${SAMPLE}_align.log                  |
					grep "overall alignment rate$"              | awk '{print $1}')
			B_CONC_PER=$(echo  ${B_CONC}"/"${PASSED_FILTERS}"*100"  | bc -l)%
			B_MULTI_PER=$(echo ${B_MULTI}"/"${PASSED_FILTERS}"*100" | bc -l)%
			B_UNAL_PER=$(echo  ${B_UNAL}"/"${PASSED_FILTERS}"*100"  | bc -l)%
			UNIQ_MAPPED=$(cat logs/deDup/${SAMPLE}_deDup.log            |
					grep "Input Reads:"                         | awk '{print $7}')
			UNIQ_MAPPED_DEDUP=$(cat logs/deDup/${SAMPLE}_deDup.log      |
					grep "Number of reads out:"                 | awk '{print $8}')
			PER_DUPS=$(echo "(1-"${UNIQ_MAPPED_DEDUP}"/"${UNIQ_MAPPED}")*100" | bc -l)%
			UNIQ_MAPPED_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
					grep "Input Reads:"                         | awk '{print $7}')
			UNIQ_MAPPED_DEDUP_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
					grep "Number of reads out:"                 | awk '{print $8}')
			PER_DUPS_SPIKE=$(echo "(1-"${UNIQ_MAPPED_DEDUP_SPIKE}"/"${UNIQ_MAPPED_SPIKE}")*100" | bc -l)%

			echo -e \
				${NAME}'\t'\
				${RAW_READS}'\t'\
				${TRIMMED_READS}'\t'\
				${PER_DIMER}'\t'\
				${RRNA}'\t'\
				${PER_RRNA}'\t'\
				${PASSED_FILTERS}'\t'\
				${B_CONC}'\t'\
				${B_MULTI}'\t'\
				${B_UNAL}'\t'\
				${B_OAP}'\t'\
				${B_CONC_PER}'\t'\
				${B_MULTI_PER}'\t'\
				${B_UNAL_PER}'\t'\
				${UNIQ_MAPPED}'\t'\
				${UNIQ_MAPPED_DEDUP}'\t'\
				${PER_DUPS}'\t'\
				${UNIQ_MAPPED_SPIKE}'\t'\
				${UNIQ_MAPPED_DEDUP_SPIKE}'\t'\
				${PER_DUPS_SPIKE} >> info/infoTable.tsv
		done
	fi

	module unload python3
	module load python

	# Making non-normalized bigWig files
	echo "bigWig"
	for FILE in BAMdeDuped/*.BAM
	do
		# making stranded bigWig with binSize 1
		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_fwd.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_fwd.bw}) \
				--binSize 1 \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagInclude 16
		fi

		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_rev.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_rev.bw}) \
				--binSize 1 \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagExclude 16
		fi

		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_fwd.${binSize}.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_fwd.${binSize}.bw}) \
				--binSize ${binSize} \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagInclude 16
		fi

		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_rev.${binSize}.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_rev.${binSize}.bw}) \
				--binSize ${binSize} \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagExclude 16 \
				--scaleFactor -1
		fi

		# making unstranded bigWig with binSize 1
		if [ ! -s "bigWig/$(basename ${FILE/.BAM/.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/.bw}) \
				--binSize 1 \
				--numberOfProcessors ${THREADS} \
				--Offset 1
		fi
	done

elif [[ ${PAIRED} == "Y" ]]; then

	### Trimming adapters and filtering rRNA reads
	echo "trimming adapters and filtering rDNA reads..."

	## Branches for either 3' UMI or both UMIs
	if [[ ${THREEP_UMI} == "Y" ]]; then
		# Branch for both UMIs
        	if [[ ${FIVEP_UMI} == "Y" ]]; then
			for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
			do
				if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]; then
					echo "trimming adapters and filtering rRNA reads for ${PAIR}"
					(fastp \
						-i fastq/${PAIR}_R1.fastq \
						-I fastq/${PAIR}_R2.fastq \
						--adapter_sequence ${ADAPTOR_1} \
						--adapter_sequence_r2 ${ADAPTOR_2} \
						--umi \
						--stdout \
						--umi_loc=per_read \
						--umi_len=${UMI_LEN} \
						--html logs/fastp/${PAIR}_fastp.html \
						-w $(echo ${THREADS}/3 | bc)\
						-c \
						--overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
					(bowtie2 \
						--fast-local \
						--un-conc trimmedFastq/${PAIR}.fastq \
						--interleaved - \
						-x ${RDNA} \
						--threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
					mv fastp.json logs/fastp/${PAIR}_fastp.json
				fi
			done
		# Branch for just 3' UMI
		else
			for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
			do
				if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]; then
					echo "trimming adapters and filtering rRNA reads for ${PAIR}"
					(fastp \
						-i fastq/${PAIR}_R1.fastq \
						-I fastq/${PAIR}_R2.fastq \
						--adapter_sequence ${ADAPTOR_1} \
						--adapter_sequence_r2 ${ADAPTOR_2} \
						--umi \
						--stdout \
						--umi_loc=read1 \
						--umi_len=${UMI_LEN} \
						--html logs/fastp/${PAIR}_fastp.html \
						-w $(echo ${THREADS}/3 | bc) \
						-c \
						--overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
					(bowtie2 \
						--fast-local \
						--un-conc trimmedFastq/${PAIR}.fastq \
						--interleaved - \
						-x ${RDNA} \
						--threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
					mv fastp.json logs/fastp/${PAIR}_fastp.json
				fi
			done
		fi
	# Branch for only 5' UMI or no UMIs
	else
		# Branch for only 5' UMI
		if [[ $FIVEP_UMI == "Y" ]]; then
			for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
			do
				if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]; then
					echo "trimming adapters and filtering rRNA reads for ${PAIR}"
					(fastp \
						-i fastq/${PAIR}_R1.fastq \
						-I fastq/${PAIR}_R2.fastq \
						--adapter_sequence ${ADAPTOR_1} \
						--adapter_sequence_r2 ${ADAPTOR_2} \
						--umi \
						--stdout \
						--umi_loc=read2 \
						--umi_len=${UMI_LEN} \
						--html logs/fastp/${PAIR}_fastp.html \
						-w $(echo ${THREADS}/3 | bc) \
						-c \
						--overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
					(bowtie2 \
						--fast-local \
						--un-conc trimmedFastq/${PAIR}.fastq \
						--interleaved - \
						-x ${RDNA} \
						--threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
					mv fastp.json logs/fastp/${PAIR}_fastp.json
				fi
			done
		# Branch for no UMI
		else
			for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
			do
				if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]; then
					echo "trimming adapters and filtering rRNA reads for ${PAIR}"
					(fastp \
						-i fastq/${PAIR}_R1.fastq \
						-I fastq/${PAIR}_R2.fastq \
						--adapter_sequence ${ADAPTOR_1} \
						--adapter_sequence_r2 ${ADAPTOR_2} \
						--stdout \
						--html logs/fastp/${PAIR}_fastp.html \
						-w $(echo ${THREADS}/3 | bc) \
						-c \
						--overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
					(bowtie2 \
						--fast-local \
						--un-conc trimmedFastq/${PAIR}.fastq \
						--interleaved - \
						-x ${RDNA} \
						--threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
					mv fastp.json logs/fastp/${PAIR}_fastp.json
				fi
			done
		fi
	fi

	### Cleaning up filenames in trimmedFastq (bowtie automatically names PE --un output)
	for FILE in trimmedFastq/*1.fastq
	do
		if [ ! -s ${FILE/.1.fastq/_R1.fastq} ]; then
			mv "${FILE}" ${FILE/.1.fastq/_R1.fastq}
		fi
	done

	for FILE in trimmedFastq/*2.fastq
	do 
		if [ ! -s ${FILE/.2.fastq/_R2.fastq} ]; then
			mv "${FILE}" ${FILE/.2.fastq/_R2.fastq}
		fi
	done

	### Aligning to spike in genome to get normalization factors

	for PAIR in $(ls trimmedFastq | sed 's/_R[1-2].*//' | uniq )
	do
		if [ ! -s "spikeBAM/${PAIR}_hg38.BAM" ]; then
			echo "aligning ${PAIR} to spike in genome"
			(bowtie2 \
				--local \
				--very-sensitive-local \
				--threads $(echo ${THREADS}/3*2 | bc) \
				--no-unal \
				--no-mixed \
				--no-discordant \
				-x "${GENOME_SPIKE}" \
				-1 "trimmedFastq/${PAIR}_R1.fastq" \
				-2 "trimmedFastq/${PAIR}_R2.fastq" \
				2> logs/spikeAlign/${PAIR}_spikeAlign.log) |
			samtools view -hS -f 2 -q ${MAPQ} |
			samtools view -b | 
			samtools sort -@ $(echo ${THREADS}/3 | bc) -o spikeBAM/${PAIR}.BAM
			samtools index spikeBAM/${PAIR}.BAM
			#perl -n -e 'print $_ if (/^\@/ || /'${SPIKE_PREFIX}'/ ) ' |
		fi
	done

	### Aligning to experimental genome

	for PAIR in $(ls trimmedFastq | sed 's/_R[1-2].*//' | uniq )
	do
		if [ ! -s "BAM/${PAIR}.BAM" ]; then
			echo "aligning ${PAIR} to experimental genome"
			(bowtie2 \
				--local \
				--sensitive-local \
				--threads $(echo ${THREADS}/3*2 | bc) \
				-x "${GENOME_EXP}" \
				-1 "trimmedFastq/${PAIR}_R1.fastq" \
				-2 "trimmedFastq/${PAIR}_R2.fastq" \
				2> logs/align/${PAIR}_align.log) |
			samtools view -bS -f 2 -q ${MAPQ} |
			samtools sort -@ ${THREADS} -o BAM/${PAIR}.BAM
			samtools index BAM/${PAIR}.BAM
			echo ""
		fi
	done

	### deduplicating with UMIs

	module unload python
	module load python3

	for FILE in BAM/*.BAM
	do
		if [ ! -s "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" ]; then
			(umi_tools dedup \
				-I "${FILE}" \
				--umi-separator=":" \
				--paired \
				-S "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
				)> "logs/deDup/$(basename ${FILE%.BAM}_deDup.log)" &&
			samtools index "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" 
		fi
	done

	## deduplicating with UMIs

	for FILE in spikeBAM/*.BAM
	do
		if [ ! -s "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" ]; then
			(umi_tools dedup \
				-I "${FILE}" \
				--umi-separator=":" \
				--paired \
				-S "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
				)> "logs/spikedeDup/$(basename ${FILE%.BAM}_deDup.log)" &&
			samtools index "spikeBAMdseDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" 
		fi
	done

	### Generating infoTable

	if [ ! -s info/infoTable.tsv ]; then
		touch info/infoTable.tsv
		echo -e \
			Name'\t'\
			RawReads'\t'\
			NonDimerReads'\t'\
			%dimer'\t'\
			insertSize'\t'\
			rRNAreads'\t'\
			%rRNA'\t'\
			passedFilters'\t'\
			bowtieConcordant'\t'\
			bowtieMulti'\t'\
			bowtieUnal'\t'\
			bowtieOverallMap%'\t'\
			bowtieConcordant%'\t'\
			bowtieMulti%'\t'\
			bowtieUnal%'\t'\
			uniqueMapped'\t'\
			uniqueMappedNondup'\t'\
			%PCRdups'\t'\
			uniqueMappedSpikein'\t'\
			uniqueMappedSpikeinNondup'\t'\
			spikeInPCRdups% >> info/infoTable.tsv

		for SAMPLE in $(ls BAM/*.BAM | sed 's/.BAM//' | sed 's/BAM\///' )
		do
			NAME=${SAMPLE}
			RAW_READS=$(cat logs/fastp/${SAMPLE}_fastp.log |
					grep "total reads:" | head -n 1 | awk '{print $3}')
			TRIMMED_READS=$(cat logs/fastp/${SAMPLE}_fastp.log |
					grep "total reads:" | tail -n 1 | awk '{print $3}')
			PER_DIMER=$(echo "(1-"${TRIMMED_READS}"/"${RAW_READS}")*100" | bc -l)%
			INSERT_SIZE=$(cat logs/fastp/${SAMPLE}_fastp.log |
					grep "Insert size peak" | awk '{print $8}')
			PASSED_FILTERS=$(cat logs/align/${SAMPLE}_align.log |
					grep "reads; of these:$" | awk '{print $1}')
			RRNA=$(echo ${TRIMMED_READS}"-"${PASSED_FILTERS} | bc )
			PER_RRNA=$(echo ${RRNA}"/"${RAW_READS}"*100" | bc -l)%
			B_CONC=$(cat logs/align/${SAMPLE}_align.log |
					grep "aligned concordantly exactly 1 time$" | awk '{print $1}')
			B_MULTI=$(cat logs/align/${SAMPLE}_align.log |
					grep "aligned concordantly >1 times$" | awk '{print $1}')
			B_UNAL=$(cat logs/align/${SAMPLE}_align.log |
					grep "aligned concordantly 0 times$" | awk '{print $1}')
			B_OAP=$(cat logs/align/${SAMPLE}_align.log |
					grep "overall alignment rate$" | awk '{print $1}')
			B_CONC_PER=$(echo ${B_CONC}"/"${PASSED_FILTERS}"*100" | bc -l)%
			B_MULTI_PER=$(echo ${B_MULTI}"/"${PASSED_FILTERS}"*100" | bc -l)%
			B_UNAL_PER=$(echo ${B_UNAL}"/"${PASSED_FILTERS}"*100" | bc -l)%
			UNIQ_MAPPED=$(cat logs/deDup/${SAMPLE}_deDup.log |
					grep "Input Reads:" | awk '{print $10}')
			UNIQ_MAPPED_DEDUP=$(cat logs/deDup/${SAMPLE}_deDup.log |
					grep "Number of reads out:" | awk '{print $8}')
			PER_DUPS=$(echo "(1-"${UNIQ_MAPPED_DEDUP}"/"${UNIQ_MAPPED}")*100" | bc -l)%
			UNIQ_MAPPED_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
					grep "Input Reads:" | awk '{print $10}')
			UNIQ_MAPPED_DEDUP_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
					grep "Number of reads out:" | awk '{print $8}')
			PER_DUPS_SPIKE=$(echo "(1-"${UNIQ_MAPPED_DEDUP_SPIKE}"/"${UNIQ_MAPPED_SPIKE}")*100" | bc -l)%

			echo -e \
				${NAME}'\t'\
				${RAW_READS}'\t'\
				${TRIMMED_READS}'\t'\
				${PER_DIMER}'\t'\
				${INSERT_SIZE}'\t'\
				${RRNA}'\t'\
				${PER_RRNA}'\t'\
				${PASSED_FILTERS}'\t'\
				${B_CONC}'\t'\
				${B_MULTI}'\t'\
				${B_UNAL}'\t'\
				${B_OAP}'\t'\
				${B_CONC_PER}'\t'\
				${B_MULTI_PER}'\t'\
				${B_UNAL_PER}'\t'\
				${UNIQ_MAPPED}'\t'\
				${UNIQ_MAPPED_DEDUP}'\t'\
				${PER_DUPS}'\t'\
				${UNIQ_MAPPED_SPIKE}'\t'\
				${UNIQ_MAPPED_DEDUP_SPIKE}'\t'\
				${PER_DUPS_SPIKE} >> info/infoTable.tsv
		done
	fi

	module unload python3
	module load python

	# Making non-normalized bigWig files

	for FILE in BAMdeDuped/*.BAM
	do
		# making stranded bigWig with binSize 1
		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_fwd.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_fwd.bw}) \
				--binSize 1 \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagInclude 82
		fi

		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_rev.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_rev.bw}) \
				--binSize 1 \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagInclude 98
		fi

		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_fwd.${binSize}.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_fwd.${binSize}.bw}) \
				--binSize ${binSize} \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagInclude 82
		fi

		if [ ! -s "bigWig/$(basename ${FILE/.BAM/_rev.${binSize}.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/_rev.${binSize}.bw}) \
				--binSize ${binSize} \
				--numberOfProcessors ${THREADS} \
				--Offset 1 \
				--samFlagInclude 98 \
				--scaleFactor -1
		fi

		# making unstranded bigWig with binSize 1
		if [ ! -s "bigWig/$(basename ${FILE/.BAM/.bw})" ]; then
			bamCoverage \
				--bam ${FILE} \
				--skipNonCoveredRegions \
				--outFileName bigWig/$(basename ${FILE/.BAM/.bw}) \
				--binSize 1 \
				--numberOfProcessors ${THREADS} \
				--Offset 1
		fi
	done
else
	echo "wrong library type"
fi

for FILE in fastq/*
do
	if [[ "${FILE}" == *.fastq ]]; then 
		gzip ${FILE} & 
	fi
done
wait

for FILE in trimmedFastq/*
do
	if [[ "${FILE}" == *.fastq ]]; then 
		gzip ${FILE} & 
	fi
done
wait

echo ""
echo "end of PROSeq mapping"

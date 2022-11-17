#!/bin/bash
#SBATCH -p standard -o matrix.log -t 8:00:00
#SBATCH -c 8 --mem=64G

module load deeptools

binSize=25
beforeRegionStartLength=2000
regionBodyLength=5000
afterRegionStartLength=2000
numberOfProcessors=4

bedFilep=/scratch/cshih8/references/annotation/gencode.hg38.protein_coding.uni.fwd.bed
bedFilem=/scratch/cshih8/references/annotation/gencode.hg38.protein_coding.uni.rev.bed

sampleFiles=$(<samples.txt)

mkdir -p QC/enrichment/transcripts
mkdir -p QC/enrichment/tss

for sample in ${sampleFiles[*]}
do
	forward=${sample}_R1_dedup_QC_plus.bw
	reverse=${sample}_R1_dedup_QC_minus.bw

	computeMatrix scale-regions \
		--scoreFileName ${forward} \
		--binSize ${binSize} \
		--regionsFileName ${bedFilep} \
		--beforeRegionStartLength ${beforeRegionStartLength} \
		--regionBodyLength ${regionBodyLength} \
		--afterRegionStartLength ${afterRegionStartLength} \
		--skipZeros \
		--outFileName matrix_f1.mat.gz \
		--missingDataAsZero \
		--samplesLabel ${sample} \
		--numberOfProcessors ${numberOfProcessors}

	computeMatrix scale-regions \
		--scoreFileName ${reverse} \
		--binSize ${binSize} \
		--regionsFileName ${bedFilem} \
		--beforeRegionStartLength ${beforeRegionStartLength} \
		--regionBodyLength ${regionBodyLength} \
		--afterRegionStartLength ${afterRegionStartLength} \
		--skipZeros \
		--outFileName matrix_f2.mat.gz \
		--missingDataAsZero \
		--samplesLabel ${sample} \
		--numberOfProcessors ${numberOfProcessors} \
		--scale -1

	computeMatrixOperations relabel -m matrix_f1.mat.gz -o matrix_f1a.mat.gz --groupLabels "${sample}"
	computeMatrixOperations relabel -m matrix_f2.mat.gz -o matrix_f2a.mat.gz --groupLabels "${sample}"
	computeMatrixOperations rbind -m matrix_f1a.mat.gz matrix_f2a.mat.gz -o matrix.mat.gz

	plotProfile \
		--matrixFile matrix.mat.gz \
		--colors red \
		--averageType mean \
		--outFileName QC/enrichment/transcripts/${sample}.pdf

	rm -f matrix*

	computeMatrix reference-point \
		--scoreFileName ${forward} \
		--binSize ${binSize} \
		--regionsFileName ${bedFilep} \
		--beforeRegionStartLength ${beforeRegionStartLength} \
		--afterRegionStartLength ${afterRegionStartLength} \
		--skipZeros \
		--outFileName matrix_f1.mat.gz \
		--missingDataAsZero \
		--samplesLabel ${sample} \
		--numberOfProcessors ${numberOfProcessors}

	computeMatrix reference-point \
		--scoreFileName ${reverse} \
		--binSize ${binSize} \
		--regionsFileName ${bedFilem} \
		--beforeRegionStartLength ${beforeRegionStartLength} \
		--afterRegionStartLength ${afterRegionStartLength} \
		--skipZeros \
		--outFileName matrix_f2.mat.gz \
		--missingDataAsZero \
		--samplesLabel ${sample} \
		--numberOfProcessors ${numberOfProcessors} \
		--scale -1

	computeMatrix reference-point \
		--scoreFileName ${forward} \
		--binSize ${binSize} \
		--regionsFileName ${bedFilem} \
		--beforeRegionStartLength ${beforeRegionStartLength} \
		--afterRegionStartLength ${afterRegionStartLength} \
		--skipZeros \
		--outFileName matrix_r1.mat.gz \
		--missingDataAsZero \
		--samplesLabel ${sample} \
		--numberOfProcessors ${numberOfProcessors}

	computeMatrix reference-point \
		--scoreFileName ${reverse} \
		--binSize ${binSize} \
		--regionsFileName ${bedFilep} \
		--beforeRegionStartLength ${beforeRegionStartLength} \
		--afterRegionStartLength ${afterRegionStartLength} \
		--skipZeros \
		--outFileName matrix_r2.mat.gz \
		--missingDataAsZero \
		--samplesLabel ${sample} \
		--numberOfProcessors ${numberOfProcessors} \
		--scale -1

	computeMatrixOperations relabel -m matrix_f1.mat.gz -o matrix_f1a.mat.gz --groupLabels "plus"
	computeMatrixOperations relabel -m matrix_f2.mat.gz -o matrix_f2a.mat.gz --groupLabels "plus"
	computeMatrixOperations relabel -m matrix_r1.mat.gz -o matrix_r1a.mat.gz --groupLabels "minus" 
	computeMatrixOperations relabel -m matrix_r2.mat.gz -o matrix_r2a.mat.gz --groupLabels "minus"

	computeMatrixOperations rbind -m matrix_f1a.mat.gz matrix_f2a.mat.gz -o matrix_f.mat.gz
	computeMatrixOperations rbind -m matrix_r1a.mat.gz matrix_r2a.mat.gz -o matrix_r.mat.gz
	computeMatrixOperations rbind -m matrix_f.mat.gz matrix_r.mat.gz     -o matrix.mat.gz

	plotProfile \
		--matrixFile matrix.mat.gz \
		--colors red green \
		--outFileName QC/enrichment/tss/${sample}.pdf

	rm -f matrix*
done

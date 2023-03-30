#!/bin/bash
#SBATCH -p standard -o bw_scaling.log -t 23:50:00
#SBATCH -c 8 --mem=64G

module load kentutils
module load bedops

chromsize="/scratch/cshih8/references/hg38/fasta/default/hg38.chrom.sizes"

mkdir -p bigWig_scaled

while read exp
do
	IFS=' ' read -r -a array <<< "${exp}"

	file=${array[0]}
	scaleFactor=${array[1]}

	echo "bigWig file: ${file}" 
	echo "scaling factor: ${scaleFactor}"

	bigWigToBedGraph ${file}.bw tmp.bedGraph

	awk -v SF="${scaleFactor}" '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$4*SF}' tmp.bedGraph > tmp_scale.bedGraph

	bedGraphToBigWig tmp_scale.bedGraph ${chromsize} bigWig_scaled/${file}_scaled.bw

	rm -f tmp.bedGraph
	rm -f tmp_scale.bedGraph

done < bw2scale.txt

# format of bw2scale
# {file name of bw} \t {scaling factor}

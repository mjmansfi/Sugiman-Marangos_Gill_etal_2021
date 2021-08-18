#!/bin/bash

SCRIPTDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
cd "${SCRIPTDIR}"/../

WORKDIR="${PWD}"
RUNMODE="docker"
THREADS=30

if [ "${RUNMODE}" == "docker" ];
then
	docker pull mjmansfi/dt:0.4.0
	INVOCATION="docker run --user $(id -u):$(id -g) --volume "${PWD}":/work --workdir /work -it mjmansfi/dt:0.4.0"
	declare -A FASTA_ARRAY
	FASTA_ARRAY["FullLength_15seqs"]="/work/data/FullSeqs_15seq.fasta"
	FASTA_ARRAY["FullLength_17seqs"]="/work/data/FullSeqs_17seq.fasta"
	FASTA_ARRAY["JustTdomains_15seq"]="/work/data/JustTdomains_15seq.fasta"
	FASTA_ARRAY["JustTdomains_17seq"]="/work/data/JustTdomains_17seq.fasta"
fi

if [ "${RUNMODE}" == "singularity" ];
then
	cd lib
	singularity pull docker://mjmansfi/dt:0.4.0 
	cd ../
	INVOCATION="singularity exec ./lib/dt_0.4.0.sif"
	declare -A FASTA_ARRAY
	FASTA_ARRAY["FullLength_15seqs"]="./data/FullSeqs_15seq.fasta"
	FASTA_ARRAY["FullLength_17seqs"]="./data/FullSeqs_17seq.fasta"
	FASTA_ARRAY["JustTdomains_15seq"]="./data/JustTdomains_15seq.fasta"
	FASTA_ARRAY["JustTdomains_17seq"]="./data/JustTdomains_17seq.fasta"

fi


mkdir -p "${PWD}"/output
for FASTA in "${!FASTA_ARRAY[@]}";
do
	mkdir -p ./output/"${FASTA}"
	for MSAALGO in mafft prank clustalo muscle;
	do
		mkdir -p ./output/"${FASTA}"/alignment
		ALIGNMENT=./output/"${FASTA}"/alignment/"${MSAALGO}".mfa
		
		# Make alignments.
		if [ "${MSAALGO}" == "mafft" ];
		then
			eval "${INVOCATION}" linsi --thread "${THREADS}" --quiet "${FASTA_ARRAY[${FASTA}]}" > "${ALIGNMENT}"
		elif [ "${MSAALGO}" == "prank" ];
		then
			eval "${INVOCATION}" prank -F -protein -o="${ALIGNMENT}" -d="${FASTA_ARRAY[${FASTA}]}"
			mv "${ALIGNMENT}".best.fas "${ALIGNMENT}"
		elif [ "${MSAALGO}" == "clustalo" ];
		then
			eval "${INVOCATION}" clustalo --threads="${THREADS}" --auto --outfile="${ALIGNMENT}" --infile="${FASTA_ARRAY[${FASTA}]}" --seqtype=Protein
		elif [ "${MSAALGO}" == "muscle" ];
		then
			eval "${INVOCATION}" muscle -in "${FASTA_ARRAY[${FASTA}]}" -out "${ALIGNMENT}"
		fi
		
		# Make trees.
		mkdir -p ./output/"${FASTA}"/raxml/"${MSAALGO}"
		eval "${INVOCATION}" raxmlHPC-PTHREADS-AVX2 -T "${THREADS}" -m PROTGAMMAAUTO -x 1992 -p 1992 -N autoMRE -f a -s "${ALIGNMENT}" -n "${MSAALGO}"
		eval "${INVOCATION}" raxmlHPC-PTHREADS-AVX2 -T "${THREADS}" -f S -s "${ALIGNMENT}" -n "${MSAALGO}".sspb -m PROTGAMMAAUTO -t RAxML_bestTree."${MSAALGO}"
		grep '^Final ML Optimization Likelihood' RAxML_info.* | awk '{print $NF}' > RAxML_likelihood."${MSAALGO}".txt
		echo "Done!"+
		eval "${INVOCATION}" Rscript ./scripts/plot_tree.R --treeFile=./output/"${FASTA}"/raxml/"${MSAALGO}"/RAxML_bipartitions."${MSAALGO}" --root=midpoint --outFile="${FASTA}"."${MSAALGO}".raxml
		mv RAxML* "${FASTA}"."${MSAALGO}".raxml* ./output/"${FASTA}"/raxml/"${MSAALGO}"

		
		mkdir -p ./output/"${FASTA}"/mrbayes/"${MSAALGO}"
		NEXUS="./output/"${FASTA}"/mrbayes/"${MSAALGO}"/"${MSAALGO}".nexus"
		eval "${INVOCATION}" python ./scripts/biopython_aligned_fasta_to_nexus.py -f "${ALIGNMENT}" -n "${NEXUS}" -m protein
		MRBAYESOPTS="begin mrbayes;\n\tset autoclose=yes nowarn=yes;\n\texecute "${NEXUS}";\n\tlset nst=6 rates=gamma;\n\tprset aamodelpr=mixed;\n\tmcmc nruns=3 ngen=1000000 samplefreq=1000 mcmcdiagn=yes diagnfreq=5000 diagnstat=avgstddev minpartfreq=0.10 relburnin=yes burninfrac=0.25 savetrees=no savebrlens=yes stoprule=yes stopval=0.01 file="${NEXUS}".mrbayes;\n\tsump;\n\tconformat=simple;\n\tsumt;\nend;\n"
		printf "${MRBAYESOPTS}" > "${NEXUS}".commands
		eval "${INVOCATION}" mb "${NEXUS}".commands > "${NEXUS}".mb.log
		eval "${INVOCATION}" Rscript ./scripts/plot_tree.R --treefile=./output/"${FASTA}"/mrbayes/"${MSAALGO}".nexus.mrbayes.con.tre --root=midpoint --treeFormat=nexus --outFile="${FASTA}"."${MSAALGO}".mrbayes
		mv "${FASTA}"."${MSAALGO}".mrbayes ./output/"${FASTA}"/mrbayes/
	done
done

#!/bin/bash

SCRIPTDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
cd "${SCRIPTDIR}"/../

WORKDIR="${PWD}"
RUNMODE="singularity"
THREADS=30

if [ "${RUNMODE}" == "docker" ];
then
	docker pull mjmansfi/dt:0.4.0
	INVOCATION="docker run --user $(id -u):$(id -g) --volume "${PWD}":/work --workdir /work -it mjmansfi/dt:0.4.0"
	declare -A FASTA_ARRAY
	FASTA_ARRAY["5seq.Full"]="/work/data/FullSeqs_5seq.fasta"
	FASTA_ARRAY["5seq.C"]="/work/data/FullSeqs_5seq.C.fasta"
	FASTA_ARRAY["5seq.T"]="/work/data/FullSeqs_5seq.T.fasta"
	FASTA_ARRAY["5seq.R"]="/work/data/FullSeqs_5seq.R.fasta"
fi

if [ "${RUNMODE}" == "singularity" ];
then
	cd lib
	singularity pull docker://mjmansfi/dt:0.4.0
	cd ../
	INVOCATION="singularity exec "${WORKDIR}"/lib/dt_0.4.0.sif"
	declare -A FASTA_ARRAY
	FASTA_ARRAY["5seq.Full"]="./data/FullSeqs_5seq.fasta"
	FASTA_ARRAY["5seq.C"]="./data/FullSeqs_5seq.C.fasta"
	FASTA_ARRAY["5seq.T"]="./data/FullSeqs_5seq.T.fasta"
	FASTA_ARRAY["5seq.R"]="./data/FullSeqs_5seq.R.fasta"
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
			eval "${INVOCATION}" linsi --thread "${THREADS}" --quiet "${FASTA_ARRAY[${FASTA}]}" > "${ALIGNMENT}" 2> /dev/null
		elif [ "${MSAALGO}" == "prank" ];
		then
			eval "${INVOCATION}" prank -F -protein -o="${ALIGNMENT}" -d="${FASTA_ARRAY[${FASTA}]}" 2> /dev/null
			mv "${ALIGNMENT}".best.fas "${ALIGNMENT}"
		elif [ "${MSAALGO}" == "clustalo" ];
		then
			eval "${INVOCATION}" clustalo --threads="${THREADS}" --auto --outfile="${ALIGNMENT}" --infile="${FASTA_ARRAY[${FASTA}]}" --seqtype=Protein 2> /dev/null
		elif [ "${MSAALGO}" == "muscle" ];
		then
			eval "${INVOCATION}" muscle -in "${FASTA_ARRAY[${FASTA}]}" -out "${ALIGNMENT}" 2> /dev/null
		fi

		# Make trees.
		mkdir -p ./output/"${FASTA}"/raxml/"${MSAALGO}"
        echo "Running RAxML for "${FASTA}" and "${MSAALGO}""
		eval "${INVOCATION}" raxmlHPC-PTHREADS-AVX2 -T "${THREADS}" -m PROTGAMMAAUTO -x 1992 -p 1992 -N autoMRE -f a -s "${ALIGNMENT}" -n "${MSAALGO}" > "${FASTA}"."${MSAALGO}".raxml.log
		eval "${INVOCATION}" raxmlHPC-PTHREADS-AVX2 -T "${THREADS}" -f S -s "${ALIGNMENT}" -n "${MSAALGO}".sspb -m PROTGAMMAAUTO -t RAxML_bestTree."${MSAALGO}" > "${FASTA}"."${MSAALGO}".raxml.sspb.log
		grep '^Final ML Optimization Likelihood' RAxML_info.* | awk '{print $NF}' > RAxML_likelihood."${MSAALGO}".txt
		mv RAxML* "${FASTA}"."${MSAALGO}".raxml* ./output/"${FASTA}"/raxml/"${MSAALGO}"


		mkdir -p ./output/"${FASTA}"/mrbayes/"${MSAALGO}"
		NEXUS="./output/"${FASTA}"/mrbayes/"${MSAALGO}"/"${MSAALGO}".nexus"
		# Note that MrBayes doesn't like fancy characters in the FASTA headers.
		cat "${ALIGNMENT}" | cut -d "/" -f 1 > tmp.fa
		eval "${INVOCATION}" python ./scripts/biopython_aligned_fasta_to_nexus.py -f tmp.fa -n "${NEXUS}" -m protein
		rm tmp.fa
        echo "Running MrBayes for "${FASTA}" and "${MSAALGO}""
		MRBAYESOPTS="begin mrbayes;\n\tset autoclose=yes nowarn=yes;\n\texecute "${NEXUS}";\n\tlset nst=6 rates=gamma;\n\tprset aamodelpr=mixed;\n\tmcmc nruns=3 ngen=1000000 samplefreq=1000 mcmcdiagn=yes diagnfreq=5000 diagnstat=avgstddev minpartfreq=0.10 relburnin=yes burninfrac=0.25 savetrees=no savebrlens=yes stoprule=yes stopval=0.01 file="${NEXUS}".mrbayes;\n\tsump;\n\tsumt burninfrac=0.25 Nruns=3 Conformat=FigTree;\nend;\n"
		printf "${MRBAYESOPTS}" > "${NEXUS}".commands
		eval "${INVOCATION}" mb "${NEXUS}".commands > "${NEXUS}".mb.log
		
		# Plot the trees. Note that directory mounting complications makes it
		# easier to keep everything in one directory, so that's why it's done this way here.
		mkdir -p ./output/summary
		ln -srf ./output/"${FASTA}"/raxml/"${MSAALGO}"/RAxML_bipartitions."${MSAALGO}" ./output/summary/"${FASTA}"."${MSAALGO}".raxml.newick
		ln -srf ./output/"${FASTA}"/mrbayes/"${MSAALGO}"/"${MSAALGO}".nexus.mrbayes.con.tre ./output/summary/"${FASTA}"."${MSAALGO}".mb.nexus
		ln -srf ./scripts/plot_tree.R ./output/summary/plot_tree.R
		cd ./output/summary
		eval "${INVOCATION}" Rscript ./plot_tree.R --treeFile="${FASTA}"."${MSAALGO}".raxml.newick --root=midpoint --treeFormat=newick --outFile="${FASTA}"."${MSAALGO}".raxml
		eval "${INVOCATION}" Rscript ./plot_tree.R --treeFile="${FASTA}"."${MSAALGO}".mb.nexus --root=midpoint --treeFormat=nexus --outFile="${FASTA}"."${MSAALGO}".mrbayes
		cd "${WORKDIR}"
		echo "Done "${FASTA}" and "${MSAALGO}"!"
	done
done

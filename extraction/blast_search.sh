#!/usr/bin/env bash

spacer_dir="$1"
ref_fa="$2" # e.g., ref/bl21de3_pcas12laci_psbk23min.fa
bin_path="$3"
err_comment="failed to find ncbi-blast bin path; please specify a proper bin path"

if [ -z "${ref_fa}" ]; then
	# If not all non-optional args given, output usage and exit
	echo "Usage:"
	echo "./blast_search.sh [spacer_dir] [reference_fasta] [bin_path (optional)]"
	exit 1
fi

if [ -z "${bin_path}" ]; then
	# Look for proper bin path if none given
	echo "Usage:"
	echo "./build_blast_db.sh [bin_path (optional)]"
	echo ""
	echo "No path to binaries given; looking along PATH variable"
	which -s blastn
	if [ $? -eq 0 ]; then
		use_cmd="blastn"
	else # if not found
		echo "${err_comment}"
		exit 1
	fi
else
	if [ -d "${bin_path}" ]; then
		use_cmd="${bin_path}/blastn"
	else
		echo "${err_comment}"
		exit 1
	fi
fi

#for all files, perform blast search
find ${spacer_dir}/*_out.fa -type f | while read nm; do
	newfn="${nm%_*}_blast_out.txt"
	echo $newfn
	${use_cmd} -db "${ref_fa}" -query "${nm}" -task "blastn" \
		-num_threads 4 -out "${newfn}" -evalue 0.0001 -outfmt 10
done


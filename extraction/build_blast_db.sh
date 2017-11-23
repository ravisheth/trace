#!/usr/bin/env bash

# Builds a BLAST database for aligning spacers

bin_path="$1"
err_comment="failed to find ncbi-blast bin path; please specify a proper bin path"

if [ -z "${bin_path}" ]; then
	# Look for proper bin path if none given
	echo "Usage:"
	echo "./build_blast_db.sh [bin_path (optional)]"
	echo ""
	echo "No path to binaries given; looking along PATH variable"
	which -s makeblastdb
	if [ $? -eq 0 ]; then
		use_cmd="makeblastdb"
	else # if not found
		echo "${err_comment}"
		exit 1
	fi
else
	if [ -d "${bin_path}" ]; then
		use_cmd="${bin_path}/makeblastdb"
	else
		echo "${err_comment}"
		exit 1
	fi
fi

# Make BLAST dbs for any .fa files in ref/ folder
find "ref/" -mindepth 1 -type f \( -iname "*.fa" -not -iname ".*" \) |
while read -r fa_fn; do
	${use_cmd} -in "${fa_fn}" -parse_seqids -dbtype nucl 
done

echo "Done"

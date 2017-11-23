#!/usr/bin/env bash

# Quick script to find raw miseq data, gunzip it to the ./data folder

search_dir="$1"
out_dir="$2"

if [ -z "${search_dir}" ]; then # if no args, output usgae and exit
	echo "Usage:"
	echo "./find_raw_miseq_data.sh [dir_to_search] [out_dir (optional)]"
	exit 1
fi

echo "Keep compressed data (y/n)? " # User input
read keep_data

if [ -z "${out_dir}" ]; then # if no out dir, then save in the dir being searched
	out_dir="${search_dir}"
fi

echo "Unzipping fastq files to ${out_dir}"
if [ "${keep_data}" == "y" ]; then
	# Keeping the .gz files
	find "${search_dir}" -mindepth 1 -type f \( -not -iname ".*" \) | \
	while read -r fn; do
		bn=`basename ${fn}`
		gzip -dvc "${fn}" > "${out_dir}/${bn%.*}"
	done
else
	# Removing the .gz files
	find ${search_dir} -mindepth 1 -type f \( -not -iname ".*" \) | \
	while read -r fn; do
		mv "${fn}" "${out_dir}"
		new_fn="${out_dir}/$( basename ${fn} )"
		gzip -dv "${new_fn}"

		enclosing_dir=`find $( dirname ${fn} ) -maxdepth 0 -empty`
		if [ ! -z "${enclosing_dir}"  ]; then # Remove old enclosing dir if empty now
			echo "Removing dir ${enclosing_dir}"
			rmdir "${enclosing_dir}"
		fi
	done
fi

echo "Done"

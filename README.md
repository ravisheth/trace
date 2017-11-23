# trace: temporal recording in arrays by crispr expansion

code from the paper **"Multiplex recording of cellular events over time on CRISPR biological tape"** Sheth RU, Yim SS, Wu FL, Wang HH. Science, 2017 (First Release)

the full paper and supplemental information can be accessed [here](http://science.sciencemag.org/lookup/doi/10.1126/science.aao0958)

## data avaliability

raw sequencing data can be found at NCBI SRA under PRJNA417866. plasmids have been deposited to Addgene; plasmid maps can be found in the _plasmid_maps_ folder.

## dependencies

* python 2.7, ipython/jupyter 5.1.0
	- pandas
	- numpy
	- matplotlib
	- seaborn
	- scipy
	- sklearn
	- _NB: Most of the above libraries are bundled together in the [Anaconda distribution](https://www.continuum.io/downloads)_
* [NCBI BLAST 2.6.0](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

we have tested this code on Mac OS X v10.11.3.

## spacer extraction and alignment workflow

unzip gzipped fastq files from a given location to a new directory within the repo (e.g., *data/my_result_dir*)
```
$ mkdir data/my_result_dir
$ ./extraction/find_unzip_raw_data.sh [dir_to_search] [out_dir (optional)]
```

make data output directory and extract spacers from the raw read data
```
$ mkdir data/my_output_dir
$ ./extraction/spacer_extraction.py [fastq_directory] [out_directory] [DR_sequence (optional)]
```

create BLAST database for determining spacer origins. this only needs to be done once. note that if the ncbi-BLAST bin is already on your path, the script can be executed without the path argument. *we only provide the reference for the main pRec/pTrig recording strain in the _ref_ folder to save space in the repo, but references for the other recording strains can be easily recreated using plasmid sequences from the _plasmid_maps_ folder.*
```
$ ./extraction/build_blast_db.sh [bin_path (optional)]
```

search the spacers against the BLAST database. again, note that if the ncbi-BLAST bin is already on your path, the script can be executed without the path argument.
```
$ ./extraction/blast_search.sh [spacer_dir] [reference_fasta] [bin_path (optional)]
```

determine unique spacers from the BLAST search results
```
$ ./extraction/unique_spacers.py [working_directory]
```

## data analysis

to get you started tinkering with the TRACE system, we have provided some example analysis code to investigate the resulting data. check out the demo notebook, where we analyze the 4 day temporal recording experiment from Fig. 2 and 3 in the Science manuscript: [_demo/trace_4day_analysis.ipynb_](demo/trace_4day_analysis.ipynb)
```
$ cd demo
$ jupyter notebook
```

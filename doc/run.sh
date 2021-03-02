# somatic
mkdir -p fastq/A_control
touch fastq/A_control/pass.txt
touch fastq/A_tumor/pass.txt

snakemake -s ./snakefile_somatic.txt --forceall --dag | dot -Tpng > dag_somatic.png

# germline
mkdir -p fastq/A_tumor
touch fastq/A_tumor/pass.txt

snakemake -s ./snakefile_germline.txt --forceall --dag | dot -Tpng > dag_germline.png


# germline
mkdir -p fastq/sampleR
touch fastq/sampleR/pass.txt

snakemake -s ./snakefile_rna.txt --forceall --dag | dot -Tpng > dag_rna.png



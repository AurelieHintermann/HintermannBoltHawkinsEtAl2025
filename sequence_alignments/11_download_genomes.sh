
download_genome() {
	local input_genome="$1"
	local output_name="$2"

	curl --output ${input_genome}.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${input_genome}/download?include_annotation_type=GENOME_FASTA" && \
	unzip -o ${input_genome}.zip
	cp ncbi_dataset/data/${input_genome}/*.fna ${output_name}_${input_genome}.fna
	rm -r ncbi_dataset
	rm ${input_genome}.zip

}

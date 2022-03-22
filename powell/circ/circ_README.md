
rule circexplorer2_annotate_prep:
	priority: 1
	input:
		GTF = annotation_gtf
	output:
		GENE_PRED = os.path.join(circexplorer2_out_dir, "mRatBN7.2_annotation.txt")
	shell:
		"""
		source ~/miniconda3/etc/profile.d/conda.sh && conda deactivate && conda activate circexplorer2 &&
		gtfToGenePred {input.GTF} temp_mRatBN7.2_annotation.txt
		cut -f1 temp_mRatBN7.2_annotation.txt > temp_firstcol.txt
		paste -d$'\\t' temp_firstcol.txt temp_mRatBN7.2_annotation.txt > {output.GENE_PRED}
		rm temp*
		"""
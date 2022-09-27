# source: https://github.com/DNAstack/plenary-resources-2020/blob/13dce9cc0131b7924a41ae80bea4b93353cae6fa/workflows/gwas.wdl

version 1.0

workflow gwas {
	input {
		File vcf
		File metadata_csv
		Int chromnum
	}

	call parse_metadata {
		input:
			metadata_csv = metadata_csv
	}

	call run_gwas {
		input:
			vcf = vcf,
			covariates = parse_metadata.covariates,
			phenotypes = parse_metadata.phenotypes,
			sex = parse_metadata.sex,
			ids = parse_metadata.ids
	}

	call create_plot {
		input:
			assoc_file = run_gwas.logistic,
			chromnum = chromnum
	}

	output {
		File logistic = run_gwas.logistic
		File manhattan_plot = create_plot.manhattan_plot
	}
}

task parse_metadata {
	input {
		File metadata_csv
	}

	command {
		parse_metadata.sh \
			-c ~{metadata_csv}
	}

	output {
		File covariates = "covariates.txt"
		File phenotypes = "phenotypes.txt"
		File sex = "sex.txt"
		File ids = "ids.txt"
	}

	runtime {
		docker: "dnastack/plink:1.9"
		cpu: 1
		memory: "3.75 GB"
		disks: "local-disk 20 HDD"
	}
}

task run_gwas {
	input {
		File vcf
		File covariates
		File phenotypes
		File sex
		File ids
	}

	# Currently query params are not stripped off when files are localized, leading to some awkward file
	# Names. The sub() is essentialyl stripping off any query params
	String output_basename = basename(sub(vcf,"\\?.*",""), ".vcf.gz")

	command {
		plink \
			--vcf ~{vcf} \
			--maf 0.10 \
			--update-ids ~{ids} \
			--make-bed \
			--out ~{output_basename}

		plink \
			--bfile ~{output_basename} \
			--update-sex ~{sex} \
			--pheno ~{phenotypes} \
			--make-bed \
			--out ~{output_basename}

		# Recode covariates to binary
		plink \
			--bfile ~{output_basename} \
			--covar ~{covariates} \
			--dummy-coding \
			--write-covar

		plink \
			--bfile ~{output_basename} \
			--logistic \
			--covar plink.cov \
			--out ~{output_basename}

		sed -i -e 's/\s\+/,/g' -e 's/^,//g' -e 's/,$//g' ~{output_basename}.assoc.logistic
	}

	output {
		File logistic = "~{output_basename}.assoc.logistic"
	}

	runtime {
		docker: "dnastack/plink:1.9"
		cpu: 4
		memory: "16 GB"
		disks: "local-disk 350 HDD"
		preemptible: 2
	}
}

task create_plot {
	input {
		File assoc_file
                Int chromnum
	}

	String assoc_basename = basename(assoc_file, ".assoc.logistic")

	command {
		manhattan_plot.py \
			-i ~{assoc_file} \
			-c ~{chromnum} \
			-o ~{assoc_basename}.png
	}

	output {
		File manhattan_plot = "~{assoc_basename}.png"
	}

	runtime {
		docker: "dnastack/plink:1.9"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk 20 HDD"
	}
}

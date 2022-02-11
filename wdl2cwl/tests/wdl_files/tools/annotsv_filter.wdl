version 1.0

task annotsvFilter {
  input {
    Boolean all_CDS = false
    File annotsv_tsv
    Float filtering_frequency = 0.05
    Boolean ignore_pass_filter = false
    String output_tsv_name = "filtered-bcftools-merged-AnnotSV.tsv"
  }

  Int space_needed_gb = 10 + round(2*size(annotsv_tsv, "GB"))
  runtime {
    memory: "4GB"
    docker: "python:3"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    python -c '
    import csv
    import sys
    input_file_name  = "~{annotsv_tsv}"
    output_file_name = "~{output_tsv_name}"
    filtering_frequency = ~{filtering_frequency}
    all_cds = ~{true="True" false="False" all_CDS}
    ignore_pass_filter = ~{true="True" false="False" ignore_pass_filter}
    with open(input_file_name, "r") as file_in, open(output_file_name, "w") as file_out:
        file_in = csv.DictReader(file_in, delimiter="\t")
        file_out = csv.DictWriter(file_out, fieldnames=file_in.fieldnames, delimiter="\t")
        file_out.writeheader()
        total_sv_count = 0
        pass_sv_count = 0
        for row in file_in:
            total_sv_count += 1
            if(row["AnnotSV type"] == "split" \
                and (row["FILTER"] == "PASS" or ignore_pass_filter) \
                and (int(row["CDS length"]) > 0 or all_cds) \
                and float(row["IMH_AF"]) < filtering_frequency
                and float(row["1000g_max_AF"]) < filtering_frequency
                and not(float(row["DGV_LOSS_Frequency"]) > filtering_frequency and "DEL" in row["SV type"])
                and not(float(row["DGV_GAIN_Frequency"]) < filtering_frequency and ("DUP" in row["SV type"] or "INS" in row["SV type"]))
                and not(("Manta" in row["ID"] and "IMPRECISE" in row["INFO"]) or (row["QUAL"] != "." and "IMPRECISE" in row["INFO"])) ):
                file_out.writerow(row)
                pass_sv_count += 1
        print("total sv count:",total_sv_count)
        print("total sv passed count:",pass_sv_count)
    '
  >>>

  output {
    File filtered_tsv = output_tsv_name
  }
}

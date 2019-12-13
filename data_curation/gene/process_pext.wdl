######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Process gnomAD pext data

workflow process_pext {
  String output_filename
  File contiglist
  String storage_bucket

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Iterate over contigs
  scatter ( contig in contigs ) {
    call format_pext {
      input:
        contig=contig[0]
    }
  }

  # Merge & sort reformatted pext data
  call mergesort_bed as mergesort_pext {
    input:
      beds=format_pext.formatted_data,
      output_filename=output_filename,
      upload_to_bucket=storage_bucket
  }
}


# Format pext data for a single chromosome
task format_pext {
  String contig
  String rCNV_bucket

  command <<<
    gsutil -m cp \
      gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.022719.tsv.bgz \
      ./
    tabix -s 1 -b 2 -e 2 -S 1 \
      all.possible.snvs.tx_annotated.022719.tsv.bgz
    /opt/rCNV2/data_curation/gene/process_pext.py \
      --contig ${contig} \
      -o pext_data.${contig}.bed.gz \
      --bgzip \
      all.possible.snvs.tx_annotated.022719.tsv.bgz
  >>>

  output {
    File formatted_data = "pext_data.${contig}.bed.gz"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:1b3910427c948c50677f89e444910bed4c43a690503a7b9efc90bfe427ad1ba8"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 50 SSD"
    bootDiskSizeGb: "20"
  }
}


# Merge & sort an array of BED files
task mergesort_bed {
  Array[File] beds
  String output_filename
  String upload_to_bucket

  command <<<
    zcat ${sep=" " beds} \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bgzip -c \
    > ${output_filename}

    gsutil -m cp ${output_filename} ${upload_to_bucket}
  >>>

  output {
    File mergesorted_bed = "${output_filename}"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:1b3910427c948c50677f89e444910bed4c43a690503a7b9efc90bfe427ad1ba8"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 50 SSD"
    bootDiskSizeGb: "20"
  }
}
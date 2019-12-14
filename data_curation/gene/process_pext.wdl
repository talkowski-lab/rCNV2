######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Process gnomAD pext data

workflow process_pext {
  File pext_gs_path
  String output_filename
  File contiglist
  String storage_bucket

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Generate pext tabix index
  call tabix_pext {
    input:
      pextfile=pext_gs_path
  }

  # Iterate over contigs
  scatter ( contig in contigs ) {
    call format_pext {
      input:
        pextfile=pext_gs_path,
        pext_idx=tabix_pext.pext_idx,
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


# Generate tabix index for pext file
task tabix_pext {
  File pextfile

  command <<<
    tabix -s 1 -b 2 -e 2 -S 1 ${pextfile}
    pextfile_basename=$( basename ${pextfile} )
    find / -name "$pextfile_basename".tbi \
    | xargs -I {} mv {} ./
  >>>

  output {
    File pext_idx = glob("*.tbi")[0]
  }

  runtime {
    docker: "talkowski/rcnv@sha256:ee72d7b02283be11db9f8642b75a71accc749de19344b067a7ffbf124cee16e5"
    preemptible: 1
  }
}


# Format pext data for a single chromosome
task format_pext {
  File pextfile
  File pext_idx
  String contig

  command <<<
    /opt/rCNV2/data_curation/gene/process_pext.py \
      --contig ${contig} \
      --pan-tissue \
      -o pext_data.${contig}.bed.gz \
      --bgzip \
      ${pextfile}
  >>>

  output {
    File formatted_data = "pext_data.${contig}.bed.gz"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:ee72d7b02283be11db9f8642b75a71accc749de19344b067a7ffbf124cee16e5"
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
    | sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 \
    | bgzip -c \
    > ${output_filename}
    tabix -f ${output_filename}

    gsutil -m cp ${output_filename}* ${upload_to_bucket}
  >>>

  output {
    File mergesorted_bed = "${output_filename}"
    File mergesorted_bed_idx = "${output_filename}.tbi"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:ee72d7b02283be11db9f8642b75a71accc749de19344b067a7ffbf124cee16e5"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 50 SSD"
    bootDiskSizeGb: "20"
  }
}


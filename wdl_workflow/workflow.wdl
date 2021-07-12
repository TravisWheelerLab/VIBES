version 1.0

workflow Pseudomonas {
  input {
    File phage_file
    Array[File] bacteria_files
  }

  call hmmbuild { input: seq_file=phage_file }
  scatter (bacteria in bacteria_files) {
    call create_table { input: hmm_file=hmmbuild.hmm_file, hmm_f_file=hmmbuild.hmm_f_file, hmm_i_file=hmmbuild.hmm_i_file, hmm_m_file=hmmbuild.hmm_m_file, hmm_p_file=hmmbuild.hmm_p_file, genome_file=bacteria }
    call format_table { input: dfam_file=create_table.scanned_dfam_file, genome_file=bacteria }
  }

  output {
    Array[File] insertion_files = format_table.insertion_file
    Array[File] occurrence_files = format_table.occurrence_file
  }
}

task hmmbuild {
  input {
    File seq_file
    String runner_name = "traviswheelerlab/pseudomonas_pipeline_runner"
    String runner_version = "latest"
  }

  output {
    File hmm_file = "${seq_file}.hmm"
    File hmm_f_file = "${seq_file}.hmm.h3f"
    File hmm_i_file = "${seq_file}.hmm.h3i"
    File hmm_m_file = "${seq_file}.hmm.h3m"
    File hmm_p_file = "${seq_file}.hmm.h3p"
  }

  command {
    perl /programs/perl/table_gen/hmmbuild_mult_seq.pl \
        --input "${seq_file}" \
        --output "${seq_file}.hmm" \
        --hmmpress \
        --cpu 11
  }

  runtime {
    docker: "${runner_name}:${runner_version}"
  }
}

task create_table {
  input {
    File hmm_file
    File hmm_f_file
    File hmm_i_file
    File hmm_m_file
    File hmm_p_file
    File genome_file
    String runner_name = "traviswheelerlab/pseudomonas_pipeline_runner"
    String runner_version = "latest"
  }

  output {
    File scanned_dfam_file = "${genome_file}.scanned.dfam"
  }

  command {
    cp "${genome_file}" "${genome_file}.dfam"
    perl /programs/perl/table_gen/dfam_tableizer.pl \
        --hmm_db "${hmm_file}" \
        --genome "${genome_file}" \
        --dfam "${genome_file}.dfam" \
        --scanned_dfam "${genome_file}.scanned.dfam"
  }

  runtime {
    docker: "${runner_name}:${runner_version}"
  }
}

task format_table {
  input {
    File dfam_file
    File genome_file
    String runner_name = "traviswheelerlab/pseudomonas_pipeline_runner"
    String runner_version = "latest"
  }

  output {
    File insertion_file = "${genome_file}.tsv"
    File occurrence_file = "${genome_file}.counts"
  }

  command {
    cp "${genome_file}" "${genome_file}.tsv"
    cp "${genome_file}" "${genome_file}.counts"
  }

  runtime {
    docker: "${runner_name}:${runner_version}"
  }
}


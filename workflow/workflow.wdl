version 1.0

workflow Pseudomonas {
  input {
    File phage_file
    Array[File] bacteria_files
  }

  call hmmbuild { input: seq_file=phage_file }
  scatter (bacteria in bacteria_files) {
    call tableize { input: hmm_file=hmmbuild.hmm_file, genome_file=bacteria }
    call format_table { input: dfam_file=tableize.dfam_file, genome_file=bacteria }
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
  }

  command {
    perl /programs/perl/table_gen/hmmbuild_mult_seq.pl \
        --input "${seq_file}" --output "${seq_file}.hmm"
  }

  runtime {
    docker: "${runner_name}:${runner_version}"
  }
}

task tableize {
  input {
    File hmm_file
    File genome_file
    String runner_name = "traviswheelerlab/pseudomonas_pipeline_runner"
    String runner_version = "latest"
  }

  output {
    File dfam_file = "${genome_file}.hmm"
  }

  command {
    cp "${genome_file}" "${genome_file}.hmm"
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


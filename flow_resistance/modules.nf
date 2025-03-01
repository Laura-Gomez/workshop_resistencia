// Antibiotic resistence

process resfinder {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_res}:/resfinder_db"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple path("resfinder_results/${sample}/pheno_table.txt"), path("resfinder_results/${sample}/ResFinder_results.txt"), path("resfinder_results/${sample}/ResFinder_results_tab.txt"),     emit: resfinder_out
  val(sample), emit: resfinder_sample_out

  script:
  """
  mkdir -p resfinder_results/${sample}

  resfinder \
	--db_path_res /resfinder_db \
	-o resfinder_results/${sample} \
	-l 0.6 -t 0.8 --acquired \
	-ifa ${fasta_genome}
 """
}


// Antibiotic resistence

process resfinderfq {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_res}:/resfinder_db"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_reads)

  output:
  tuple path("resfinderfq_results/${sample}/pheno_table.txt"), path("resfinderfq_results/${sample}/ResFinder_results.txt"), path("resfinderfq_results/${sample}/ResFinder_results_tab.txt"),     emit: resfinder_out
  val(sample), emit: resfinder_sample_out

  script:
  """
  mkdir -p resfinderfq_results/${sample}

  resfinder \
        --db_path_res /resfinder_db \
        -o resfinderfq_results/${sample} \
        -l 0.6 -t 0.8 --acquired \
        -ifq ${fastp_reads[0]} ${fastp_reads[1]}
 """
}


// Run Virulence Finder

process virulencefinder {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_vir}:/virfinder_db"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("virulence/${sample}/results.txt"),
		path("virulence/${sample}/results_tab.tsv"),    emit: virfinder_out


  script:
  """
  mkdir -p virulence/${sample}

  virulencefinder.py \
	-i ${fasta_genome} \
	-p /virfinder_db \
	-o virulence/${sample} \
	--mincov 0.6 -t 0.9 \
	-x

 """

}


// Run Virulence Finder FROM FASTQ

process virulencefinderfq {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.db_vir}:/virfinder_db"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_reads)


  output:
  tuple val(sample), path("virulence_fastq/${sample}/results.txt"),
                path("virulence_fastq/${sample}/results_tab.tsv"),    emit: virfinder_out


  script:
  """
  mkdir -p virulence_fastq/${sample}

  virulencefinder.py \
        -i ${fastp_reads[0]} ${fastp_reads[1]} \
        -p /virfinder_db \
        -o virulence_fastq/${sample} \
        --mincov 0.6 -t 0.9 \
        -x

 """

}


// Resistance Gene Identifier (RGI) software 

process rgi {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)


  output:
  tuple val(sample), path("virulence_rgi/${sample}/results.txt"), emit: rgi_out


  script:
  """
  mkdir -p virulence_rgi/${sample}

  rgi load --card_json ${params.card_db} --local

  rgi main \
	--input_sequence ${fasta_genome}  \
	--output_file virulence_rgi/${sample}/results \
	--local \
	--low_quality
 """

}


// SCCMEC
process sccmec {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("sccmec/${sample}/sccmec.regions.blastn.tsv"),
		path("sccmec/${sample}/sccmec.regions.details.tsv"), 
		path("sccmec/${sample}/sccmec.targets.blastn.tsv"),
		path("sccmec/${sample}/sccmec.targets.details.tsv"),
		path("sccmec/${sample}/sccmec.tsv"), emit: sccmec_out


  script:
  """
  mkdir -p sccmec/${sample}
  sccmec --input ${fasta_genome} --outdir sccmec/${sample}/
  
  """

}



// GENE ANNOTATION
process prokka {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fasta_genome)

  output:
  tuple val(sample), path("prokka/${sample}/genome_annotation.faa"), 
		path("prokka/${sample}/genome_annotation.log"), emit: prokaa_faa

  script:
  """
  prokka --outdir prokka/${sample} --prefix genome_annotation ${fasta_genome} 
  """

}

// HIDDEN MARKOV MODELS SCAN
process hmmscan {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-resistance:public'
  containerOptions "-v ${params.resfam_dir}:/resfam"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(prokka_faa), path(prokka_log)

  output:
  tuple val(sample), path("hmmscan/${sample}_table.txt"), emit: hmmscan_out

  script:
  """
  mkdir -p hmmscan/
  hmmscan -o hmmscan/${sample}.report --cpu ${params.ncrs} --tblout hmmscan/${sample}_table.txt /resfam/${params.resfam_name} ${prokka_faa}

  """
}





process validInput {
  cache 'lenient'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("${sample}_R1.fastq.gz"),\
         path("${sample}_R2.fastq.gz"),  emit: valid_out, optional: true

  script:
  """
  if gzip -t ${reads[0]}; then
     if gzip -t ${reads[1]}; then
	  cp ${reads[0]} ${sample}_R1.fastq.gz
	  cp ${reads[1]} ${sample}_R2.fastq.gz
     fi
  fi
  """
}


// Remove adapters
process fastp {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(reads0), path(reads1)

  output:
  tuple val(sample), path("fastp/${sample}_R1.fastq.gz"),\
	 path("fastp/${sample}_R2.fastq.gz"),  emit: fastp_out
	 path("fastp/${sample}.html"), emit: fastp_html
         path("fastp/${sample}.json"), emit: fastp_json

  script:
  """
    mkdir -p fastp/
    fastp \
	--in1 ${reads0} \
	--in2 ${reads1} \
	--out1 fastp/${sample}_R1.fastq.gz \
	--out2 fastp/${sample}_R2.fastq.gz \
	--html fastp/${sample}.html \
  	--json fastp/${sample}.json
  """
}

// Continue only if high quality reads exist

process filterQual {
  cache 'lenient'

  input:
  tuple val(sample), path(fastp_1), path(fastp_2)

  output:
  tuple val(sample), path("fastp_filtered/${sample}_R1.fastq.gz"),\
         path("fastp_filtered/${sample}_R2.fastq.gz"),  emit: fastp_filtered_out, optional: true

  script:
  """
  mkdir -p fastp_filtered

  gunzip -c ${fastp_1} > ${sample}_R1.fastq

  SIZE=\$(du ${sample}_R1.fastq | cut -f1)
  INT_SIZE=\$(printf "%.0f" \"\${SIZE}\")

  if (( \"\${INT_SIZE}\" > 0)) ; then
    cp ${fastp_1} fastp_filtered/${sample}_R1.fastq.gz
    cp ${fastp_2} fastp_filtered/${sample}_R2.fastq.gz
  fi
  """
}

// Align against genome

process bwa_1 {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  containerOptions "-v ${params.refdir_8325}:/ref"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("bwa/${sample}_1.bam"), 
		path("bwa/${sample}_1.bam.bai"), emit: bwa_out

  script:
  """

  mkdir -p bwa

  bwa mem  /ref/${params.refname_8325} \
        ${fastp_data_R1} \
        ${fastp_data_R2} \
	-t 12 |
  samtools sort -o bwa/${sample}_1.bam
  samtools index bwa/${sample}_1.bam bwa/${sample}_1.bam.bai
 """
}

process bwa_2 {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  containerOptions "-v ${params.refdir_N315}:/ref"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("bwa/${sample}_2.bam"),
		path("bwa/${sample}_2.bam.bai"), emit: bwa_out

  script:
  """

  mkdir -p bwa

  bwa mem  /ref/${params.refname_N315} \
        ${fastp_data_R1} \
        ${fastp_data_R2} \
        -t 12 |
  samtools sort -o bwa/${sample}_2.bam
  samtools index bwa/${sample}_2.bam bwa/${sample}_2.bam.bai
 """
}

process bwa_3 {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  containerOptions "-v ${params.refdir_CA12}:/ref"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("bwa/${sample}_3.bam"), 
		path("bwa/${sample}_3.bam.bai"), emit: bwa_out

  script:
  """

  mkdir -p bwa

  bwa mem  /ref/${params.refname_CA12} \
        ${fastp_data_R1} \
        ${fastp_data_R2} \
        -t 12 |
  samtools sort -o bwa/${sample}_3.bam
  samtools index bwa/${sample}_3.bam bwa/${sample}_3.bam.bai
 """
}


// Coverage statistics

process mosdepth {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(bwa), path(bwa_index)

  output:
  tuple val(sample), \
        path("mosdepth/${bwa}.mosdepth.global.dist.txt"), emit: mosdepth_out

  script:
  """

  mkdir -p mosdepth

  mosdepth  mosdepth/${bwa} ${bwa}
 """
}


// Assemble

process spades {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("spades/${sample}_scaffolds.fasta"), emit: spades_out

  script:
  """

  mkdir -p spades
  
  spades.py \
	--isolate \
        -1 ${fastp_data_R1} \
        -2 ${fastp_data_R2} \
        -o spades

 cp spades/scaffolds.fasta spades/${sample}_scaffolds.fasta

 """
}


// Assemble comparison
process quast_1 {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  containerOptions "-v ${params.refdir_8325}:/ref"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold)

  output:
  tuple path("quast/${sample}_8325/report.html"),  path("quast/${sample}_8325/report.tsv"),   emit: align_out

  script:
  """
  mkdir -p quast/${sample}_8325
  
  FILE2=\$(echo ${spades_scaffold} | sed s/.fasta/_8325.fasta/)
  cp -r ${spades_scaffold} \${FILE2}

  quast.py \${FILE2} \
	-r /ref/${params.refname_8325} \
	-g /ref/${params.gff_8325} \
	-o quast/${sample}_8325
  """
}


// Assemble comparison
process quast_2 {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  containerOptions "-v ${params.refdir_N315}:/ref"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold)

  output:
  tuple path("quast/${sample}_N315/report.html"),  path("quast/${sample}_N315/report.tsv"),   emit: align_out

  script:
  """
  mkdir -p quast/${sample}_N315

  FILE2=\$(echo ${spades_scaffold} | sed s/.fasta/_N315.fasta/)
  cp -r ${spades_scaffold} \${FILE2}
  quast.py \${FILE2} \
        -r /ref/${params.refname_N315} \
        -g /ref/${params.gff_N315} \
        -o quast/${sample}_N315
  """
}


// Assemble comparison
process quast_3 {
  cache 'lenient'
  container 'laugoro/workshop-inmegen-assembly:public'
  containerOptions "-v ${params.refdir_CA12}:/ref"
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold)

  output:
  tuple path("quast/${sample}_CA12/report.html"),  path("quast/${sample}_CA12/report.tsv"),   emit: align_out

  script:
  """
  mkdir -p quast/${sample}_CA12

  FILE2=\$(echo ${spades_scaffold} | sed s/.fasta/_CA12.fasta/)
  cp -r ${spades_scaffold} \${FILE2}
  quast.py \${FILE2} \
        -r /ref/${params.refname_CA12} \
        -g /ref/${params.gff_CA12} \
        -o quast/${sample}_CA12
  """
}

// Quality analysis FASTQ

process fastqc {
   cache 'lenient'
   container 'laugoro/workshop-inmegen-assembly:public'
   publishDir params.out, mode:'copy'

   input:
   tuple val(sample_id), path(reads0), path(reads1)

   output:
   path("fastqc/*"),     emit: fastqc_out

   script:
   """
   mkdir -p fastqc
   fastqc -o fastqc ${reads0} ${reads1}
   
   """
}


// Multiqc: fastqc, fastp, quast
process multiqc {
   cache 'lenient'
   publishDir params.out, mode:'copy'

   input:
   val(sample_id)
   path(dir_all)

   output:
   path("multiqc/*"),     emit: multiqc_out

   script:
   """
   mkdir -p multiqc
   multiqc -o multiqc/ ${dir_all}

   """
}


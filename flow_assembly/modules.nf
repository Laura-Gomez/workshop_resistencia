

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
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("alignment/${sample}_1.bam"), 
		path("alignment/${sample}_1.bam.bai"), emit: bwa_out

  script:
  """

  mkdir -p alignment

  bwa mem  ${params.ref_8325} \
        ${fastp_data_R1} \
        ${fastp_data_R2} \
	-t 12 |
  samtools sort -o alignment/${sample}_1.bam
  samtools index alignment/${sample}_1.bam alignment/${sample}_1.bam.bai
 """
}

process bwa_2 {
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("alignment/${sample}_2.bam"),
		path("alignment/${sample}_2.bam.bai"), emit: bwa_out

  script:
  """

  mkdir -p alignment

  bwa mem  ${params.ref_N315} \
        ${fastp_data_R1} \
        ${fastp_data_R2} \
        -t 12 |
  samtools sort -o alignment/${sample}_2.bam
  samtools index alignment/${sample}_2.bam alignment/${sample}_2.bam.bai
 """
}

process bwa_3 {
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("alignment/${sample}_3.bam"), 
		path("alignment/${sample}_3.bam.bai"), emit: bwa_out

  script:
  """

  mkdir -p alignment

  bwa mem  ${params.ref_CA12} \
        ${fastp_data_R1} \
        ${fastp_data_R2} \
        -t 12 |
  samtools sort -o alignment/${sample}_3.bam
  samtools index alignment/${sample}_3.bam alignment/${sample}_3.bam.bai
 """
}


// Coverage statistics

process mosdepth {
  cache 'lenient'
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
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("spades_scaffolds/${sample}_scaffolds.fasta"),\
         path("spades_scaffolds/${sample}_bacterial_scaffolds.fasta"),	emit: spades_out

  script:
  """

  mkdir -p spades
  mkdir -p spades_bacterial
  mkdir -p spades_scaffolds

 spades.py \
	-1 ${fastp_data_R1} \
	-2 ${fastp_data_R2} \
	-o spades
  cp spades/scaffolds.fasta spades_scaffolds/${sample}_scaffolds.fasta

  spades.py \
	--isolate \
        -1 ${fastp_data_R1} \
        -2 ${fastp_data_R2} \
        -o spades_bacterial
  cp spades_bacterial/scaffolds.fasta spades_scaffolds/${sample}_bacterial_scaffolds.fasta

 """
}


// Assemble comparison
process quast_1 {
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)

  output:
  tuple path("quast/${sample}_8325/report.html"),  path("quast/${sample}_8325/report.tsv"),   emit: align_out

  script:
  """
  mkdir -p quast/${sample}
  
  FILE1=\$(echo ${spades_scaffold_bacterial} | sed s/fasta/_8325.fasta/)
  FILE2=\$(echo ${spades_scaffold} | sed s/fasta/_8325.fasta/)
  cp ${spades_scaffold_bacterial} \${FILE1}
  cp ${spades_scaffold} \${FILE2}

  quast.py \${FILE1} \${FILE2} \
	-r ${params.ref_8325} \
	-g ${params.gff_8325} \
	-o quast/${sample}_8325
  """
}


// Assemble comparison
process quast_2 {
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)

  output:
  tuple path("quast/${sample}_N315/report.html"),  path("quast/${sample}_N315/report.tsv"),   emit: align_out

  script:
  """
  mkdir -p quast/${sample}

  FILE1=\$(echo ${spades_scaffold_bacterial} | sed s/fasta/_N315.fasta/)
  FILE2=\$(echo ${spades_scaffold} | sed s/fasta/_N315.fasta/)
  cp ${spades_scaffold_bacterial} \${FILE1}
  cp ${spades_scaffold} \${FILE2}
  quast.py \${FILE1} \${FILE2} \
        -r ${params.ref_N315} \
        -g ${params.gff_N315} \
        -o quast/${sample}_N315
  """
}


// Assemble comparison
process quast_3 {
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)

  output:
  tuple path("quast/${sample}_CA12/report.html"),  path("quast/${sample}_CA12/report.tsv"),   emit: align_out

  script:
  """
  mkdir -p quast/${sample}

  FILE1=\$(echo ${spades_scaffold_bacterial} | sed s/fasta/_CA12.fasta/)
  FILE2=\$(echo ${spades_scaffold} | sed s/fasta/_CA12.fasta/)
  cp ${spades_scaffold_bacterial} \${FILE1}
  cp ${spades_scaffold} \${FILE2}
  quast.py \${FILE1} \${FILE2} \
        -r ${params.ref_CA12} \
        -g ${params.gff_CA12} \
        -o quast/${sample}_CA12
  """
}



// Antibiotic resistence

process resfinder {
  conda '/scratch/home/lgomez/miniconda3/envs/resfinder'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)

  output:
  tuple path("resfinder_results/${sample}/pheno_table.txt"), path("resfinder_results/${sample}/ResFinder_results.txt"), path("resfinder_results/${sample}/ResFinder_results_tab.txt"),     emit: resfinder_out
  tuple path("resfinder_results_bacterial/${sample}/pheno_table.txt"), path("resfinder_results_bacterial/${sample}/ResFinder_results.txt"), path("resfinder_results_bacterial/${sample}/ResFinder_results_tab.txt"),     emit: resfinder_bacterial_out
  val(sample), emit: resfinder_sample_out

  script:
  """
  mkdir -p resfinder_results/${sample}
  mkdir -p resfinder_results_bacterial/${sample}

  python -m resfinder \
	--db_path_res ${params.db_res} \
	-o resfinder_results/${sample} \
	-l 0.6 -t 0.8 --acquired \
	-ifa ${spades_scaffold}

 python -m resfinder \
	--db_path_res ${params.db_res} \
        -o resfinder_results_bacterial/${sample} \
        -l 0.6 -t 0.8 --acquired \
        -ifa ${spades_scaffold_bacterial}
 """
}


// Antibiotic resistence

process resfinderfq {
  conda '/scratch/home/lgomez/miniconda3/envs/resfinder'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple path("resfinderfq_results/${sample}/pheno_table.txt"), path("resfinderfq_results/${sample}/ResFinder_results.txt"), path("resfinderfq_results/${sample}/ResFinder_results_tab.txt"),     emit: resfinder_out
  val(sample), emit: resfinder_sample_out

  script:
  """
  mkdir -p resfinderfq_results/${sample}

  python -m resfinder \
        --db_path_res ${params.db_res} \
        -o resfinderfq_results/${sample} \
        -l 0.6 -t 0.8 --acquired \
        -ifq ${fastp_data_R1} ${fastp_data_R2}
 """
}


// Run Virulence Finder

process virulencefinder {
  conda '/scratch/home/lgomez/miniconda3/envs/virfind'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)

  output:
  tuple val(sample), path("virulence/${sample}/results.txt"),
		path("virulence/${sample}/results_tab.tsv"),    emit: virfinder_out


  script:
  """
  mkdir -p virulence/${sample}

  virulencefinder.py \
	-i ${spades_scaffold_bacterial} \
	-p ${params.db_vir} \
	-o virulence/${sample} \
	--mincov 0.6 -t 0.9 \
	-x

 """

}


// Run Virulence Finder FROM FASTQ

process virulencefinderfq {
  conda '/scratch/home/lgomez/miniconda3/envs/virfind'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)


  output:
  tuple val(sample), path("virulence_fastq/${sample}/results.txt"),
                path("virulence_fastq/${sample}/results_tab.tsv"),    emit: virfinder_out


  script:
  """
  mkdir -p virulence_fastq/${sample}

  virulencefinder.py \
        -i ${fastp_data_R1} ${fastp_data_R2} \
        -p ${params.db_vir} \
        -o virulence_fastq/${sample} \
        --mincov 0.6 -t 0.9 \
        -x

 """

}


// Resistance Gene Identifier (RGI) software 

process rgi {
  conda '/scratch/home/lgomez/miniconda3/envs/rgi'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)


  output:
  tuple val(sample), path("virulence_rgi/${sample}/results.txt"), emit: rgi_out


  script:
  """
  mkdir -p virulence_rgi/${sample}

  rgi load --card_json ${params.card_db} --local

  rgi main \
	--input_sequence ${spades_scaffold_bacterial}  \
	--output_file virulence_rgi/${sample}/results \
	--local \
	--low_quality
 """

}


// SCCMEC
process sccmec {
  conda '/scratch/home/lgomez/miniconda3/envs/sccmec'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)

  output:
  tuple val(sample), path("sccmec/${sample}/sccmec.regions.blastn.tsv"),
		path("sccmec/${sample}/sccmec.regions.details.tsv"), 
		path("sccmec/${sample}/sccmec.targets.blastn.tsv"),
		path("sccmec/${sample}/sccmec.targets.details.tsv"),
		path("sccmec/${sample}/sccmec.tsv"), emit: sccmec_out


  script:
  """
  mkdir -p sccmec/${sample}
  sccmec --input ${spades_scaffold} --outdir sccmec/${sample}/
  
  """

}



// GENE ANNOTATION
process prokka {
conda '/scratch/home/lgomez/miniconda3/envs/resfam'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(spades_scaffold_bacterial), path(spades_scaffold)

  output:
  tuple val(sample), path("prokka/${sample}/genome_annotation.faa"), 
		path("prokka/${sample}/genome_annotation.log"), emit: prokaa_faa

  script:
  """
  prokka --outdir prokka/${sample} --prefix genome_annotation ${spades_scaffold} 
  """

}

// HIDDEN MARKOV MODELS SCAN
process hmmscan {
  conda '/scratch/home/lgomez/miniconda3/envs/resfam'
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(prokka_faa), path(prokka_log)

  output:
  tuple val(sample), path("hmmscan/${sample}_table.txt"), emit: hmmscan_out

  script:
  """
  mkdir -p hmmscan/
  hmmscan -o hmmscan/${sample}.report --cpu ${params.ncrs} --tblout hmmscan/${sample}_table.txt ${params.resfam} ${prokka_faa}

  """
}

// META-MARC
process mmarc {
  cache 'lenient'
  publishDir params.out, mode:'copy'
  
  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)
  
  output:
  tuple val(sample), path("mmarc/${sample}.tblout.scan"), emit: mmarc_out
  
  script:
  """
  mkdir -p mmarc/
  /scratch/home/lgomez/meta-marc/bin/mmarc \
	--i1 ${fastp_data_R1} \
	--i2 ${fastp_data_R2} \
	-o mmarc -f ${sample} \
	-d -l 3 -m -t ${params.ncrs}
 
  """
}


// AMR-plus-plus
process amrpp {
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(fastp_data_R1), path(fastp_data_R2)

  output:
  tuple val(sample), path("mmarc/${sample}.tblout.scan"), emit: mmarc_out

  script:
  """
  nextflow run main_AMR++.nf -profile conda

  mkdir -p mmarc/

  /scratch/home/lgomez/meta-marc/bin/mmarc \
        --i1 ${fastp_data_R1} \
        --i2 ${fastp_data_R2} \
        -o mmarc -f ${sample} \
        -d -l 3 -m -t ${params.ncrs}

  """
}



// Quality analysis FASTQ

process fastqc {
   cache 'lenient'
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


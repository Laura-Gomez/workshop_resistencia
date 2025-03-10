
// Remove adapters
process fastp {
  cache 'lenient'
  publishDir params.out, mode:'copy'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("fastp/${sample}_R1.fastq.gz"),\
         path("fastp/${sample}_R2.fastq.gz"),  emit: fastp_out
         path("fastp/${sample}.html"), emit: fastp_html
         path("fastp/${sample}.json"), emit: fastp_json

  script:
  """
    mkdir -p fastp/
    fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 fastp/${sample}_R1.fastq.gz \
        --out2 fastp/${sample}_R2.fastq.gz \
        --html fastp/${sample}.html \
        --json fastp/${sample}.json
  """
}

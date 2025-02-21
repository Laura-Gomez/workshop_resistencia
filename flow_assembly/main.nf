#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validInput;
	  fastqc as fastqcPrev;
	  fastqc as fastqcFinal;
	  fastp;
	  filterQual;
	  bwa_1;
	  bwa_2;
	  bwa_3;
	  mosdepth as mosdepth_1;
	  mosdepth as mosdepth_2;
	  mosdepth as mosdepth_3;
	  spades;
	  quast_1;
	  quast_2;
	  quast_3;
	  multiqc
 	} from "./modules.nf"	


workflow {
	
    data_fq = Channel.fromFilePairs("${params.reads}")
             .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }


   // Remove adapters
   validInput(data_fq)
   fastqcPrev(validInput.out.valid_out)
   fastp(validInput.out.valid_out)
   filterQual(fastp.out.fastp_out)
   bwa_1(filterQual.out.fastp_filtered_out)
   bwa_2(filterQual.out.fastp_filtered_out)
   bwa_3(filterQual.out.fastp_filtered_out)
   mosdepth_1(bwa_1.out.bwa_out)
   mosdepth_2(bwa_2.out.bwa_out)
   mosdepth_3(bwa_3.out.bwa_out)
   spades(filterQual.out.fastp_filtered_out)
   quast_1(spades.out.spades_out)
   quast_2(spades.out.spades_out)
   quast_3(spades.out.spades_out)
   multiqc(spades.out.spades_out.collect(),"${params.out}")
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastp;
	  bwa;
	} from "./modules.nf"	


workflow {
	
    data_fq = Channel.fromFilePairs("${params.reads}")
             .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }


   // Remove adapters
   fastp(data_fq)

   // Align against reference genome
   bwa(fastp.out.fastp_out)
 }

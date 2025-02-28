#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { resfinder;
	  resfinderfq;
	  virulencefinder;
	  virulencefinderfq;
	  rgi;
	  sccmec;
	  prokka;
	  hmmscan;
 	} from "./modules.nf"	


workflow {
	
  fastas_dataset = Channel
                .fromPath(params.fastas)
                .map { file -> tuple(file.baseName, file) }

   fastqs_dataset = Channel.fromFilePairs("${params.reads}")
             .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }

   resfinder(fastas_dataset)
   resfinderfq(fastqs_dataset)
   virulencefinder(fastas_dataset)
   virulencefinderfq(fastqs_dataset)
//   rgi(fastas_dataset)
   sccmec(fastas_dataset)
   prokka(fastas_dataset)
   hmmscan(prokka.out.prokaa_faa)
}

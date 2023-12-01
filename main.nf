#!/usr/bin/env nextflow
import nextflow.splitter.CsvSplitter

def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { spyGenomics } from  './modules/StrepGenomics.nf'


//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------


if(!params.readFoler) {
    throw new Exception("Missing parameter params.build")
  }

if(!params.reference) {
    throw new Exception("Missing parameter params.annotationName")
  }
if(!params.results) {
    throw new Exception("Missing parameter params.results")
  }

if (params.local){
    /*sample_ch =   Channel
      .fromPath([params.reads + '/*.fastq', params.reads + '/*.fastq.gz', params.reads + '/*.fq.gz'])
      .splitFastq( by : params.splitChunk, file:true  )
      */
      reads_ch = Channel.fromFilePairs([params.readFoler + '/*_{1,2}.fastq', params.readFoler + '/*_{1,2}.fastq.gz', params.readFoler + '/*_{1,2}.fq.gz', params.readFoler + '/*_R{1,2}_001.fastq.gz'])

      reads_qc = Channel.fromPath("${params.readFoler}/*", checkIfExists: true) 
} else {
    input = fetchRunAccessions(params.sraAccession)
    sample_ch = Channel.fromList(input)
}
//--------------------------------------
// Process the workflow
//-------------------------------------

workflow {
    spyGenomics(reads_ch,reads_qc)
}

#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process fastqQulity {
    
    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'

    publishDir "${params.result}/QC", mode: 'copy'   
 
    input:
    path reads

    output:
    path("${sample_id}_fastqc_out")

    script:
    sample_id = reads.getSimpleName()
    template 'fastqc.bash'
}

process multiqc {


    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'
    publishDir "${params.result}/multiqc", mode: 'copy'
    input:
    path(qc_files)

    output:
    path("*multiqc*")

    script:
    """
    multiqc ${qc_files}
    """

}

process trimming {
    

    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val("${sample_id}"), path("${sample_id}_paired*.fq.gz"), emit: trimmed_fastqs


    script:
    template 'trimming.bash'

}

process spades_assembly {


    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'
    publishDir "${params.result}/Fasta", mode: 'copy', pattern: "*.fasta"
    publishDir "${params.result}/Assembly", mode: 'copy', pattern: "${sample_id}"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}"), emit: assemblies
    tuple val("${sample_id}"), path("${sample_id}.fasta"), emit: fasta
    path("${sample_id}.fasta"), emit: mlst_fasta
    val(sample_id), emit: sample_id

    script:
    template 'spades.bash'
    
}


process quast {
    

    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'

    input:
    path(fasta)
    path(reference)
    val(sample_id)

    output:
    path("${sample_id}")

    script:
    """
    quast.py ${fasta} -r  ${reference} -o ${sample_id}
    """
}

process quastMultiqc {

    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'
    publishDir "${params.result}/Qusts_multiqc", mode: 'copy'
    input:
    path(qc_files)

    output:
    path("*multiqc*")

    script:
    """
    multiqc ${qc_files}
    """

}
process mlst_check {

    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'
    publishDir "${params.result}/MLST", mode: 'copy', pattern: "*.txt"
    
    input:
    val(scheme)
    path(fasta_files)

    output:
    path("MLST.txt")

    script:

    template 'mlst.bash'

}

process abricate {


    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'
    publishDir "${params.result}/AMR_VF", mode: 'copy', pattern: "*.txt"

    input:
    path(fasta)

    output:
    path("*txt*")

    script:

    template 'abricate.bash'

}

process prokka {

    publishDir "${params.result}/Prokka", mode: 'copy', pattern: "${sample_id}"
    conda 'home2/sybah/.conda/envs/Enterobacteriaceae_analysis'
    publishDir "${params.result}/GFF", mode: 'copy', pattern: "*.gff"

    input:
    tuple val(sample_id), path(assembly)
    val(genus)
    val(spps)

    output:
    path("${sample_id}"), emit: annotation
    path("${sample_id}.gff"), emit: gff

    script:

    template 'prokka.bash'
    
}



process emmTyping {
   
   publishDir "${params.result}/EMM", mode: 'copy'
   input:
   tuple val(sample_id), path(assembly)
   path(emmDb)

   output:
   path("${sample_id}")

   script:
   template 'emmTyping.bash'

}



workflow spyGenomics {
    
    take:
    reads_ch
    reads_qc

    main:
    
    qc = fastqQulity(reads_qc)

    multiqc(qc.collect())

    trimmed_reads = trimming(reads_ch)

    spades = spades_assembly(trimmed_reads.trimmed_fastqs)

    quast_out = quast(spades.mlst_fasta, params.reference, spades.sample_id)

    quastMultiqc(quast_out.collect())

    prokka_out = prokka(spades.fasta, 'Streptococcus', 'pyogenes')

    mlst_result = mlst_check('spyogenes',spades.mlst_fasta.collectFile())
    
    abricate_out = abricate(spades.mlst_fasta.collect())

    emm_typing = emmTyping(spades.fasta, params.emmdb)
}




#!/usr/bin/env nextflow
nextflow.enable.dsl=2


reads_ch = Channel.fromFilePairs([params.readFoler + '/*_{1,2}.fastq', params.readFoler + '/*_{1,2}.fastq.gz', params.readFoler + '/*_{1,2}.fq.gz'])

reads_qc = Channel.fromPath("${params.readFoler}/*", checkIfExists: true) 

params.result = 'Results'

// fastQC quality control process

process fastqQulity {

    publishDir "${params.result}/QC", mode: 'copy'
    container = 'biocontainers/fastqc:v0.11.9_cv7'
    
    input:
    path reads

    output:
    path("${sample_id}_fastqc_out")

    script:
    sample_id = reads.getSimpleName()
    template 'fastqc.bash'
}

process multiqc {

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
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val("${sample_id}"), path("${sample_id}_paired*.fq.gz"), emit: trimmed_fastqs


    script:
    template 'trimming.bash'

}

process spades_assembly {

    //container = 'staphb/spades'

    publishDir "${params.result}/Fasta", mode: 'copy', pattern: "*.fasta"
    publishDir "${params.result}/Assembly", mode: 'copy', pattern: "${sample_id}"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}"), emit: assemblies
    tuple val("${sample_id}"), path("${sample_id}_contigs.fasta"), emit: fasta
    path("${sample_id}_contigs.fasta"), emit: mlst_fasta

    script:

    """
    spades.py -k 21,33,55,77 --pe1-1 ${reads[0]}  --pe1-2 ${reads[1]} --careful -o ${sample_id}
    cp ${sample_id}/contigs.fasta ${sample_id}_contigs.fasta
    """
}

process mlst_check {

    publishDir "${params.result}/MLST", mode: 'copy', pattern: "*.txt"
    container = 'staphb/mlst'
    
    input:
    val(scheme)
    path(fasta_files)

    output:
    path("MLST.txt")

    script:

    """
    mlst --scheme ${scheme} --novel Novel_alleles --nopath  ${fasta_files} > MLST.txt
    """


}

process abricate {

    container = 'staphb/abricate'

    input:
    path(fasta)

    output:

    script:

    """
    abricate --db vfdb --quiet ${fasta} > AMR.txt
    abricate --summary AMR.txt > summary.txt
    """
}

process prokka {

    publishDir "${params.result}/Prokka", mode: 'copy', pattern: "${sample_id}"

    container = 'staphb/prokka'

    input:
    tuple val(sample_id), path(assembly)
    val(genus)
    val(spps)

    output:
    path("${sample_id}"), emit: annotation
    path("${sample_id}/${sample_id}.gff"), emit: gff

    script:

    """
    prokka --outdir ${sample_id} --addgenes  --genus ${genus} --species ${spps}  --usegenus --force  --prefix ${sample_id} ${assembly}
    
    """
}


workflow {
    qc = fastqQulity(reads_qc)

    multiqc(qc.collect())

    trimmed_reads = trimming(reads_ch)

    spades = spades_assembly(trimmed_reads.trimmed_fastqs)
    spades.mlst_fasta.view()
    prokka_out = prokka(spades.fasta, 'Streptococcus', 'pyogenes')

    mlst_result = mlst_check('spyogenes',spades.mlst_fasta.collect())
    
    abricate(spades.mlst_fasta.collect())

}




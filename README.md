# Streptococcus genomics workflow

This workflow is for the analysis for Streptococcus from genome assembly, annotation, mlst and associated analysis.
It is develop using the nextflow pipeline using docker containers that contain all the required software tools.

***<p align=center>Nextflow workflow</p>*** 

```mermaid
flowchart TD
    p0((Channel.fromFilePairs))
    p1((Channel.fromPath))
    p2[spyGenomics:fastqQulity]
    p3([collect])
    p4[spyGenomics:multiqc]
    p5(( ))
    p6[spyGenomics:trimming]
    p7[spyGenomics:spades_assembly]
    p8(( ))
    p9(( ))
    p10[spyGenomics:quast]
    p11([collect])
    p12[spyGenomics:quastMultiqc]
    p13(( ))
    p14(( ))
    p15(( ))
    p16[spyGenomics:prokka]
    p17(( ))
    p18(( ))
    p19([collectFile])
    p20(( ))
    p21[spyGenomics:mlst_check]
    p22(( ))
    p23([collect])
    p24[spyGenomics:abricate]
    p25(( ))
    p26(( ))
    p27[spyGenomics:emmTyping]
    p0 -->|reads_ch| p6
    p1 -->|reads_qc| p2
    p2 --> p3
    p3 --> p4
    p4 --> p5
    p6 --> p7
    p7 --> p8
    p7 --> p16
    p7 --> p10
    p7 -->|sample_id| p10
    p9 -->|reference| p10
    p10 --> p11
    p11 --> p12
    p12 --> p13
    p14 -->|genus| p16
    p15 -->|spps| p16
    p16 --> p18
    p16 --> p17
    p7 --> p19
    p19 --> p21
    p20 -->|scheme| p21
    p21 --> p22
    p7 --> p23
    p23 --> p24
    p24 --> p25
    p7 --> p27
    p26 -->|emmDb| p27
```

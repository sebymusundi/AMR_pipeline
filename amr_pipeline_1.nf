#!/usr/bin/env nextflow
// Define parameters 


params.reads="/home/bioinfo/bioinfo/rachel_data/final_data/*.fastq"
params.outdir="/home/bioinfo/bioinfo/rachel_data/results"

println "reads: $params.reads"
println "outdir: $params.outdir"

log.info """\
			  AMR PIPELINE
			------------------------------------------
			reads:	 "$params.reads"
			outdir: "$params.outdir"
			------------------------------------------
			"""
		.stripIndent()

// quality assessment using fastqc 

process fastqc{
		container 'biocontainers/fastqc:v0.11.9_cv8'
		publishDir "${params.outdir}" , mode: 'copy'
		containerOptions = "--user root"
		

		input:
		path reads
		
		output:
		path "fastqc_out"

		script:
		"""
		mkdir fastqc_out
		fastqc ${reads} -o fastqc_out --threads 8
		"""
}

// Adapter trimming using porechop

process porechop {
		container 'biocontainers/porechop:v0.2.4dfsg-1-deb_cv1 ' 
		containerOptions = "--user root"
		publishDir "${params.outdir}", mode: 'copy'
		
		
		
		input:
		path reads 

		output:
		tuple val(reads.baseName), path("${reads.baseName}_trimmed.fastq")

		script:
		"""
		porechop -i ${reads} -o ${reads.baseName}_trimmed.fastq 
		"""
}

// quality control of trimmed reads 
process fastqc_trimmed{
						container 'biocontainers/fastqc:v0.11.9_cv8'
						publishDir "${params.outdir}" , mode: 'copy'
						
						
						input:
						tuple val(x), path(porechop_trimmed)

						output:
						path 'q_trimmed'

						script:
						"""
						mkdir q_trimmed
						fastqc ${porechop_trimmed} -o q_trimmed --thread 8
						"""
}

// filter reads using filtlong 

process filtlong {
					container 'nanozoo/filtlong:latest'
					publishDir "${params.outdir}" , mode: 'copy'

					input:
					tuple val(x), path(trimmed_reads) 

					output:
					tuple val(x), path("${trimmed_reads.baseName}_filtered.fastq")

					script:
					"""
					filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000  ${trimmed_reads} > ${trimmed_reads.baseName}_filtered.fastq
					"""
}

// de novo assembly using flye 

process flye_assembly {
					 container 'nanozoo/flye:2.9.1--bba1957'
					 publishDir "${params.outdir}", mode: 'copy'
					 memory '24 GB'


					 input:
					 tuple val(x), path(assembled)

					 output:
					 tuple val(x),path("assembly/${assembled.baseName}.fasta")
					

					 script:
					"""
					 mkdir -p assembly
                     flye --nano-raw ${assembled} -g 4.4m -o assembly/${assembled.baseName}
                     mv assembly/${assembled.baseName}/assembly.fasta assembly/${assembled.baseName}.fasta

					"""

}



				
// Quality control for the assembled genome using Quast and Bandage

process quast  {
					container 'staphb/quast:latest'
					publishDir "${params.outdir}", mode: 'copy'
						
						
					input:
					tuple val(x),path(assembled_file)

					output:
					path 'quast_output'

					script:
					"""
					mkdir quast_output 
					quast.py -o quast_ouput ${assembled_file}  --thread 8
					"""
}


//create index for assembled genome 

process reads_index {
					container 'pegi3s/bwa:latest'
					containerOptions = "--user root"
					publishDir "${params.outdir}" , mode: 'copy'
					cache true


					input:
					tuple val(x),path(assembled_file)

					output:
					
					tuple val(x), path ("${assembled_file}.bwt"), 
					path ("${assembled_file}.pac"),
   					path ("${assembled_file}.ann"),
    				path ("${assembled_file}.amb"),
   					path ("${assembled_file}.sa")



					script:
					"""
					bwa index ${assembled_file} 

					"""
 }

// map reads to the indexed genome 

process align_reads{
					container 'pegi3s/bwa:latest' 
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true
					

					input:
					tuple val(x), path(filtered_reads)
					tuple val(assembled_file), path ("${assembled_file}.bwt"), path ("${assembled_file}.pac"), path ("${assembled_file}.ann"),path ("${assembled_file}.amb"), path ("${assembled_file}.sa")
					

					output:
				   	tuple val(x), path("${filtered_reads.baseName}.sam")
					
					script:
					"""	
					mkdir -p assembled
                    mv ${assembled_file}.bwt  ${assembled_file}.pac ${assembled_file}.ann ${assembled_file}.amb ${assembled_file}.sa assembled/
                    bwa mem -x ont2d assembled/${assembled_file} ${filtered_reads} > ${filtered_reads.baseName}.sam

					"""
}
 

process samtools {
					container 'staphb/samtools:latest'
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true

					input:
					tuple val(x),path(sam_file)

					output:
					tuple val(x), path("${sam_file.baseName}.sorted.bam"), path("${sam_file.baseName}.sorted.bam.bai")

					script:
					"""
					samtools view -O BAM ${sam_file} -o ${sam_file.baseName}.bam
					samtools sort ${sam_file.baseName}.bam -o ${sam_file.baseName}.sorted.bam -O BAM
					samtools index ${sam_file.baseName}.sorted.bam 
					"""
}

//  Quality control statistics for bam file using samtools stats 

process samtools_stats {
					container 'staphb/samtools:latest'
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true

					input:
					tuple val(x), path(bam_file), path(bam_file_1)

					output:
					path "${bam_file.baseName}.stats"

					script:
					"""
					samtools stats ${bam_file} > ${bam_file.baseName}.stats
					"""
}


//process racon
process racon  {
					container 'biocontainers/racon:v1.3.2-1b1-deb_cv1'
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions ="--user root"
					cache true

					input:
					tuple val(x), path(assembled), path(filtered_reads),path(sam_file)
					

					output:				
					tuple val(x), path("${filtered_reads.baseName}.racon.fasta")
					
					script:
					"""
					racon ${filtered_reads} ${sam_file} ${assembled} > ${filtered_reads.baseName}.racon.fasta
				
 					"""

}

//output:


//process medaka


process medaka  {
					container 'nanozoo/medaka:1.7.2--aa54076'
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true

					input:
					tuple val(x), path(filtered_reads), path(racon_fasta)


					output:
					tuple val(x), path("medaka_output/${racon_fasta.baseName}.medaka.fasta") 


					script:
					"""
					mkdir -p medaka_output
					medaka_consensus -m r941_min_high_g360 -i ${filtered_reads} -d ${racon_fasta} -o medaka_output -f
					mv medaka_output/consensus.fasta  medaka_output/${racon_fasta.baseName}.medaka.fasta
 					"""

}





// process circlator 
process circlator {
					container 'sangerpathogens/circlator:latest'
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"

					input:
					path assembled_file 
					path filtered_reads


					output:
					path "circlator_output"

					script:
					"""
					mkdir circlator_output
					circlator all ${assembled_file}  ${filtered_reads} circlator_output
 					"""
}


// process  prokka  genome annotation

process prokka  {
					container 'staphb/prokka:latest'
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true

					input:
					tuple val(x), path(assembled_file)


					output:
					tuple val(x), path("prokka_output/${assembled_file.baseName}.faa") 


					script:
					"""
					mkdir -p prokka_output
					prokka --kingdom Bacteria --outdir prokka_output ${assembled_file} --force
					mv medaka_output/PROKKA_03302023.faa prokka_output/${assembled_file.baseName}.faa
 					"""

}




// Determine AMR genes 
process abricate {
					container 'staphb/abricate:latest'
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"

					input:
					tuple val(x), path(assembled_file)

					output:
					tuple path("${assembled_file.baseName}.ncbiamr.csv"), 
					path("${assembled_file.baseName}.resfinderamr.csv"), 
					path("${assembled_file.baseName}.cardamr.csv"), 
					path("${assembled_file.baseName}.vfdbamr.csv"), 
					path("${assembled_file.baseName}.megaresamr.csv"),
					path("${assembled_file.baseName}.plasmidfinderamr.csv")
					
					script:
					"""
					abricate ${assembled_file} --db ncbi --csv > ${assembled_file.baseName}.ncbiamr.csv 
					abricate ${assembled_file} --db card --csv > ${assembled_file.baseName}.cardamr.csv 
					abricate ${assembled_file} --db resfinder --csv >${assembled_file.baseName}.resfinderamr.csv  
					abricate ${assembled_file} --db vfdb --csv > ${assembled_file.baseName}.vfdbamr.csv 
					abricate ${assembled_file} --db megares --csv > ${assembled_file.baseName}.megaresamr.csv 
					abricate ${assembled_file} --db plasmidfinder --csv >${assembled_file.baseName}.plasmidfinderamr.csv
					"""

}

process amr_finder {
					container 'staphb/ncbi-amrfinderplus:latest'
					publishDir "${params.outdir}", mode: 'copy'
					containerOptions = "--user root"

					input:
					tuple val(x), path(assembled_file)

					output:
					path "${assembled_file.baseName}.amr.tsv"
					
					script:
					"""
					amrfinder -n ${assembled_file} --output ${assembled_file.baseName}.amr.tsv
					"""
}


workflow{
my_reads=Channel.fromPath("$params.reads")
//my_reads.view()
fastqc_ch=fastqc(my_reads)
porechop_ch=porechop(my_reads)
fastqc_trimmed(porechop_ch)
filt_long_ch=filtlong(porechop_ch)
//filt_long_ch.view()
flye_ch=flye_assembly(filt_long_ch)
//flye_ch.view()
quast(flye_ch)
indexed_ch=reads_index(flye_ch)
//indexed_ch.view()
align_ch = align_reads(filt_long_ch, indexed_ch)
samtools_align=samtools(align_ch)
//samtools_align.view()
samtools_stats(samtools_align)
//filt_long_ch.view()
//align_ch.view()
//flye_ch.view()
flye_filt=flye_ch.combine(filt_long_ch, by: 0)
pre_racon=flye_filt.combine(align_ch, by:0)
racon_ch=racon(pre_racon)
//racon_ch.view()
//pre_medaka=filt_long_ch.combine(racon_ch, by:0)
//pre_medaka.view()
//medaka_ch=medaka(pre_medaka)
//medaka_ch.view()
//abricate(medaka_ch)
//amr_finder(medaka_ch)
//prokka(medaka_ch)
}


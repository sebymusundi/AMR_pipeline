#!/usr/bin/env nextflow
// Define parameters 


params.reads="/home/bioinfo/bioinfo/rachel_data/final_data/sample_08.fastq"
params.outdir="/home/bioinfo/bioinfo/rachel_data/results"

println "reads: $params.reads"
println "outdir: $params.outdir"


// quality assessment using fastqc 

process fastqc{
		container 'biocontainers/fastqc:v0.11.9_cv8'
		publishDir "${params.outdir}" , mode: 'copy'
		containerOptions = "--user root"
		

		input:
		path reads
		
		output:
		path("quality_control/${reads.baseName}_fastqc.html")
        path("quality_control/${reads.baseName}_fastqc.zip")

		script:
		"""
		mkdir quality_control
		fastqc ${reads} -o quality_control
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
		mkdir porechop_results
		porechop -i ${reads} -o ${reads.baseName}_trimmed.fastq 
		"""
}

// quality control of trimmed reads 
process fastqc_trimmed{
						container 'biocontainers/fastqc:v0.11.9_cv8'
						publishDir "${params.outdir}" , mode: 'copy'
						
						
						input:
						tuple val(sample_id), path(porechop_trimmed)

						output:
						tuple val(sample_id), 
                		path("quality_control_trimmed/${porechop_trimmed.baseName}_fastqc.html"), 
                		path("quality_control_trimmed/${porechop_trimmed.baseName}_fastqc.zip")

						script:
						"""
						mkdir quality_control_trimmed
						fastqc ${porechop_trimmed} -o quality_control_trimmed
						"""
}

// filter reads using filtlong 

process filtlong {
					container 'nanozoo/filtlong:latest'
					publishDir "${params.outdir}" , mode: 'copy'

					input:
					tuple val(x), path(trimmed_reads) 

					output:
					tuple val(x), path("filtlong_results/${trimmed_reads.baseName}_filtered.fastq")

					script:
					"""
					mkdir filtlong_results
					filtlong --min_length 500 --keep_percent 90 --target_bases 500000000  ${trimmed_reads} > filtlong_results/${trimmed_reads.baseName}_filtered.fastq
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



				
// Quality control for the assembled genome using Quast 

process quast  {
					container 'staphb/quast:latest'
					publishDir "${params.outdir}", mode: 'copy'
						
						
					input:
					tuple val(x),path(assembled_file)

					output:
					path("quast_output/${assembled_file.baseName}*")			

					script:
					"""
					mkdir quast_output 
					quast.py ${assembled_file}  -o quast_output/${assembled_file.baseName}
					"""
}

// Assessing the level of completeness of the genome using BUSCO

process BUSCO  {
					container 'ezlabgva/busco:v5.4.7_cv1'
					publishDir "${params.outdir}", mode: 'copy'
						
						
					input:
					tuple val(x),path(assembled_file)

					output:
					path("Busco_output/${assembled_file.baseName}*")			

					script:
					"""
					mkdir Busco_output 
					busco -m genome -i ${assembled_file} -l  proteobacteria_odb10  -o Busco_output/${assembled_file.baseName}
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
				   	tuple val(x), path("sam_files/${filtered_reads.baseName}.sam")
					
					script:
					"""	
					mkdir -p assembled
					mkdir -p sam_files
                    mv ${assembled_file}.bwt  ${assembled_file}.pac ${assembled_file}.ann ${assembled_file}.amb ${assembled_file}.sa assembled/
                    bwa mem -x ont2d assembled/${assembled_file} ${filtered_reads} > sam_files/${filtered_reads.baseName}.sam

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
					tuple val(x), path("bam_files/${sam_file.baseName}.sorted.bam"), path("bam_files/${sam_file.baseName}.sorted.bam.bai")

					script:
					"""
					mkdir -p bam_files
					samtools view -O BAM ${sam_file} -o ${sam_file.baseName}.bam
					samtools sort ${sam_file.baseName}.bam -o bam_files/${sam_file.baseName}.sorted.bam -O BAM
					samtools index bam_files/${sam_file.baseName}.sorted.bam 
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
					container 'nanozoo/racon:latest '
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions ="--user root"
					cache true
					

					input:
					tuple val(x), path(assembled), path(filtered_reads),path(sam_file)
					

					output:				
					tuple val(x), path("racon_output/${filtered_reads.baseName}.racon.fasta")
					
					script:
					"""
					mkdir racon_output
					racon ${filtered_reads} ${sam_file} ${assembled} > racon_output/${filtered_reads.baseName}.racon.fasta
				
 					"""

}

//create index for assembled genome 

process reads_index_second {
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

					script:
					"""
					mkdir index_files_2
					bwa index  ${assembled_file} index_files_2
					"""
 }

// map reads to the indexed genome 

process align_reads_second {
					container 'pegi3s/bwa:latest' 
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true
					

					input:
					tuple val(x), path(filtered_reads)
					tuple val(assembled_file), path ("${assembled_file}.bwt"), path ("${assembled_file}.pac"), path ("${assembled_file}.ann"),path ("${assembled_file}.amb"), path ("${assembled_file}.sa")
					

					output:
				   	tuple val(x), path("sam_files_2/${filtered_reads.baseName}.sam")
					
					script:
					"""	
					mkdir -p assembled_2
					mkdir -p sam_files_2
                    mv ${assembled_file}.bwt  ${assembled_file}.pac ${assembled_file}.ann ${assembled_file}.amb ${assembled_file}.sa assembled_2/
                    bwa mem -x ont2d assembled_2/${assembled_file} ${filtered_reads} > sam_files_2/${filtered_reads.baseName}.sam

					"""
}
 

process racon_second  {
					container 'nanozoo/racon:latest '
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions ="--user root"
					cache true
					

					input:
					tuple val(x), path(assembled), path(filtered_reads),path(sam_file)
					

					output:				
					tuple val(x), path("racon_output_2/${filtered_reads.baseName}.racon.fasta")
					
					script:
					"""
					mkdir racon_output_2
					racon ${filtered_reads} ${sam_file} ${assembled} > racon_output_2/${filtered_reads.baseName}.racon.fasta
				
 					"""

}

//create index for assembled genome 

process reads_index_third {
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

					script:
					"""
					mkdir index_files_3
					bwa index  ${assembled_file} index_files_3
					"""
 }

// map reads to the indexed genome 

process align_reads_third{
					container 'pegi3s/bwa:latest' 
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true
					

					input:
					tuple val(x), path(filtered_reads)
					tuple val(assembled_file), path ("${assembled_file}.bwt"), path ("${assembled_file}.pac"), path ("${assembled_file}.ann"),path ("${assembled_file}.amb"), path ("${assembled_file}.sa")
					

					output:
				   	tuple val(x), path("sam_files_3/${filtered_reads.baseName}.sam")
					
					script:
					"""	
					mkdir -p assembled_3
					mkdir -p sam_files_3
                    mv ${assembled_file}.bwt  ${assembled_file}.pac ${assembled_file}.ann ${assembled_file}.amb ${assembled_file}.sa assembled_3/
                    bwa mem -x ont2d assembled_3/${assembled_file} ${filtered_reads} > sam_files_3/${filtered_reads.baseName}.sam

					"""
}
 


process racon_third {
					container 'nanozoo/racon:latest '
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions ="--user root"
					cache true
					

					input:
					tuple val(x), path(assembled), path(filtered_reads),path(sam_file)
					

					output:				
					tuple val(x), path("racon_output_3/${filtered_reads.baseName}.racon.fasta")
					
					script:
					"""
					mkdir racon_output_3
					racon ${filtered_reads} ${sam_file} ${assembled} > racon_output_3/${filtered_reads.baseName}.racon.fasta
				
 					"""

}

//create index for assembled genome 

process reads_index_forth {
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

					script:
					"""
					mkdir index_files_4
					bwa index  ${assembled_file} index_files_4
					"""
 }

// map reads to the indexed genome 

process align_reads_forth {
					container 'pegi3s/bwa:latest' 
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions = "--user root"
					cache true
					

					input:
					tuple val(x), path(filtered_reads)
					tuple val(assembled_file), path ("${assembled_file}.bwt"), path ("${assembled_file}.pac"), path ("${assembled_file}.ann"),path ("${assembled_file}.amb"), path ("${assembled_file}.sa")
					

					output:
				   	tuple val(x), path("sam_files_4/${filtered_reads.baseName}.sam")
					
					script:
					"""	
					mkdir -p assembled_4
					mkdir -p sam_files_4
                    mv ${assembled_file}.bwt  ${assembled_file}.pac ${assembled_file}.ann ${assembled_file}.amb ${assembled_file}.sa assembled_4/
                    bwa mem -x ont2d assembled_4/${assembled_file} ${filtered_reads} > sam_files_4/${filtered_reads.baseName}.sam

					"""
}
 
process racon_forth {
					container 'nanozoo/racon:latest '
					publishDir "${params.outdir}" , mode: 'copy'
					containerOptions ="--user root"
					cache true
					

					input:
					tuple val(x), path(assembled), path(filtered_reads),path(sam_file)
					

					output:				
					tuple val(x), path("racon_output_4/${filtered_reads.baseName}.racon.fasta")
					
					script:
					"""
					mkdir racon_output_4
					racon ${filtered_reads} ${sam_file} ${assembled} > racon_output_4/${filtered_reads.baseName}.racon.fasta
				
 					"""

}
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



//process homopolish 

process homopolish {
					container 'staphb/homopolish:latest'
					publishDir "${params.outdir}", mode: 'copy'
					containerOptions = "--user root"
					cache true

					input:
					tuple val(x), path(assembled_file)

					output:
					tuple val(x), path("homopolish_output/${x}/*")

					script:
					"""
					mkdir -p homopolish_output
					homopolish polish -a ${assembled_file} -s bacteria.msh -m R9.4.pkl -o homopolish_output/${x}/				
					"""

}

// Quality control for the assembled genome using Quast 

process quast_polished{
					container 'staphb/quast:latest'
					publishDir "${params.outdir}", mode: 'copy'
						
						
					input:
					tuple val(x),path(assembled_file)

					output:
					path("quast_output_polished/${x}*")			

					script:
					"""
					mkdir quast_output_polished
					quast.py ${assembled_file}  -o quast_output_polished/${x}
					"""
}

// Assessing the level of completeness of the genome using BUSCO

process busco_polished {
					container 'ezlabgva/busco:v5.4.7_cv1'
					publishDir "${params.outdir}", mode: 'copy'
						
						
					input:
					tuple val(x),path(assembled_file)

					output:
					path("Busco_output_polished/${x}*")			

					script:
					"""
					mkdir Busco_output_polished 
					busco -m genome -i ${assembled_file} -l  proteobacteria_odb10  -o Busco_output_polished/${x}
					"""
}



// Determine individual species using MLST

process MLST {
				container 'staphb/mlst:latest'
				publishDir "${params.outdir}" , mode: 'copy'
				containerOptions = "--user root"
				cache true

				input:
				tuple val(x), path(assembled_file)
				
				output:
				tuple val(x), path("mlst_output/${x}_mlst.csv")

				script:
				"""
				mkdir mlst_output
				mlst ${assembled_file} --csv > mlst_output/${x}_mlst.csv
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
					mv medaka_output/PROKKA_*.faa prokka_output/${assembled_file.baseName}.faa
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
					tuple path("abricate_results/${assembled_file.baseName}.ncbiamr.csv"), 
					path("abricate_results/${assembled_file.baseName}.resfinderamr.csv"), 
					path("abricate_results/${assembled_file.baseName}.cardamr.csv"), 
					path("abricate_results/${assembled_file.baseName}.vfdbamr.csv"), 
					path("abricate_results/${assembled_file.baseName}.megaresamr.csv"),
					path("abricate_results/${assembled_file.baseName}.plasmidfinderamr.csv")
					
					script:
					"""
					mkdir abricate_results/
					abricate ${assembled_file} --db ncbi --csv > abricate_results/${assembled_file.baseName}.ncbiamr.csv 
					abricate ${assembled_file} --db card --csv > abricate_results/${assembled_file.baseName}.cardamr.csv 
					abricate ${assembled_file} --db resfinder --csv > abricate_results/${assembled_file.baseName}.resfinderamr.csv  
					abricate ${assembled_file} --db vfdb --csv > abricate_results/${assembled_file.baseName}.vfdbamr.csv 
					abricate ${assembled_file} --db megares --csv > abricate_results/${assembled_file.baseName}.megaresamr.csv 
					abricate ${assembled_file} --db plasmidfinder --csv >abricate_results/${assembled_file.baseName}.plasmidfinderamr.csv
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
// Retrieve nanopore reads that have been basecalled and barcoded using GUPPY

my_reads=Channel.fromPath("$params.reads")

//my_reads.view()

// Assess the quality of reads 

//fastqc_ch=fastqc(my_reads)

// Trim reads to remove adapters 
porechop_ch=porechop(my_reads)

// Check the quality of the trimmed reads 

//fastqc_trimmed(porechop_ch)

// Filter out to remove long reads 

filt_long_ch=filtlong(porechop_ch)
//filt_long_ch.view()

// Assemble the genome using flye 

flye_ch=flye_assembly(filt_long_ch)
//flye_ch.view()

// Asssess the quality of the draft assesmbly 
quast(flye_ch)
BUSCO(flye_ch)

// index the individual reads 
indexed_ch=reads_index(flye_ch)
//indexed_ch.view()

// Map the reads against the indexed draft 
align_ch = align_reads(filt_long_ch, indexed_ch)
//align_ch.view()

// Generate the BAM file
samtools_align=samtools(align_ch)
//samtools_align.view()

// Retrive the number of mapped reads 
samtools_stats(samtools_align)


// Run the first racon process

flye_filt=flye_ch.combine(filt_long_ch, by: 0)
pre_racon=flye_filt.combine(align_ch, by:0)
racon_ch=racon(pre_racon)
//racon_ch.view()


// index reads for second RACON PROCESS

index_second=reads_index_second(racon_ch)

// map reads for second racon process
align_ch_second = align_reads_second(filt_long_ch, index_second)


// The second racon process

racon_second_ch=racon_ch.combine(filt_long_ch, by: 0)
pre_second_racon_ch=racon_second_ch.combine(align_ch_second, by:0)
second_racon=racon_second(pre_second_racon_ch)
//second_racon.view()


// index reads for third RACON PROCESS 
 index_third=reads_index_third(second_racon)

// map reads for third racon process
 align_ch_third= align_reads_third(filt_long_ch, index_third)


racon_third_ch=second_racon.combine(filt_long_ch, by: 0)
pre_third_racon_ch=racon_third_ch.combine(align_ch_third, by:0)
third_racon=racon_third(pre_third_racon_ch)
//third_racon.view()



//Run the fourth RACON PROCESS
// INDEX READS
 index_forth=reads_index_forth(third_racon)

// map reads for fourth racon process
 align_ch_forth= align_reads_forth(filt_long_ch, index_forth)



racon_forth_ch=third_racon.combine(filt_long_ch, by:0)
pre_forth_racon_ch=racon_forth_ch.combine(align_ch_forth, by:0)
forth_racon=racon_forth(pre_forth_racon_ch)
//fourth_racon.view()

// Polishing the genome once using Medaka 

pre_medaka=filt_long_ch.combine(forth_racon, by:0)
//pre_medaka.view()
medaka_ch=medaka(pre_medaka)
//medaka_ch.view()


// use homopolish to ppolish the medaka genome

medaka_polished=homopolish(medaka_ch)
//medaka_polished.view()


// Provide quality metrics for the polished genome 

quast_polished(medaka_polished)

// Assess the completeness of the genome 

busco_polished(medaka_polished)

// Assign individual sequence types for individual genomes 

MLST(medaka_ch)

//Predict antimicrobial resistance genes from the provided genome 

abricate(medaka_ch)

//Predict antimicrobial resistance genes from the provided genome 

//amr_finder(medaka_ch)

//prokka(medaka_ch)
}


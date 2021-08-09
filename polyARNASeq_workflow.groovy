title: "polyA_RNA-Seq analysis pipeline"

DataDir =  "/n/scratch3/users/z/zc79/Dppa2_4/raw_data"

def branches = [

     sample1      : ["${DataDir}/sample1.R1.fastq.gz",
                     "${DataDir}/sample1.R2.fastq.gz"],

     sample2	  : ["${DataDir}/sample2.R1.fastq.gz",
                     "${DataDir}/sample2.R2.fastq.gz"]

]

Fastqc ={
      doc title: "QC using fastQC",
      desc: "Using fastQC to generate quality control reports",      
      author: "zchen0423@gmail.com"

    var SOUT : "./FastQC"
    var threds: 2
    def basename = new File(input).getName().prefix;

	exec """
      module load fastqc/0.11.5 &&
      mkdir -p ${SOUT}/${branch.name} &&
	  fastqc --threads ${threds} --extract  --outdir ${SOUT}/${branch.name} $inputs
      ""","Fastqc"
	forward inputs
}

Trimgalor = {
      def basename1 = new File(input1).getName().prefix.replaceAll(".fastq","");
      def basename2 = new File(input2).getName().prefix.replaceAll(".fastq","");

    produce("./Trimgalor/${branch.name}_trimgalor/${basename1}_val_1.fq.gz","./Trimgalor/${branch.name}_trimgalor/${basename2}_val_2.fq.gz"){
     exec """
       module load trimgalore/0.4.5 &&
       mkdir -p Trimgalor &&
        mkdir -p Trimgalor/${branch.name}_trimgalor &&
       trim_galore  --paired  --fastqc --output_dir "Trimgalor/${branch.name}_trimgalor" $inputs
     ""","Trimgalor"
    }
}

hisat2_map = {
       def hist2_index="/n/groups/zhanglab/zychen/genomes_annotations/mm10/GRC38.p6/mm10"
       def nbcpu = 5
       def outdir = "${branch.name}_hisat2"
       def alignment_splicesites = "/n/groups/zhanglab/zychen/genomes_annotations/mm10/annotations/Mus_musculus.GRCm38.97.fixed.splice.site.txt"
       branch.hisat2_dir = outdir

       produce("${outdir}/${branch.name}.bam"){
              exec """
              module load gcc boost hisat2/2.1.0 samtools/1.9 &&              
              mkdir -p ${outdir} &&
              hisat2 -x ${hist2_index} 
                     -1 ${input1} 
                     -2 ${input2} 
                     --known-splicesite-infile $alignment_splicesites 
                     --sp 1000,1000
                     --no-mixed                     
                     --no-discordant
                     -p ${nbcpu} 
                     --un-conc-gz ${outdir}/${branch.name}.unmapped.hisat2.gz
                     --met-stderr 
                     --new-summary 
                     --summary-file ${outdir}/${branch.name}.hisat2_summary.txt 
                     | samtools view -bS -F 4 -F 8 -F 256 - > ${outdir}/${branch.name}.bam 
              ""","hisat2_map"
       }
}

sortAndIndex = {
       def nbcpu = 2
       produce("${input.prefix}.sorted.bam"){
              exec """
              module load gcc boost samtools/1.9 &&
              samtools sort $input.bam -@ ${nbcpu} -o ${input.prefix}.sorted.bam &&
              samtools index ${input.prefix}.sorted.bam
              ""","sortAndIndex"
       }
}

stringTie = {
       def outdir = "stringTie/${branch.name}"
       def gtf = "/n/groups/zhanglab/zychen/genomes_annotations/mm10/annotations/mm10_plus_ercc92_fixed.gtf"
       def nbcpu = 5
       produce("$outdir/${branch.name}.gene_abund.txt"){
       exec """
                module load gcc boost stringtie/1.3.3b &&
                mkdir  -p $outdir &&
                stringtie $input 
                -o $outdir/${branch.name}_transcripts.gtf 
                -v 
                -G $gtf 
                -A $outdir/${branch.name}.gene_abund.txt
                -C $outdir/${branch.name}.cov_refs.gtf 
                -b $outdir/${branch.name}_ballgown 
                -e
                -p $nbcpu
       ""","stringTie"   
    }
    forward input
}

bamCoverage = {
   def basename = new File(input).getName().prefix
   produce("bigwigs/${basename}.bw"){
     exec """
       mkdir -p bigwigs &&
       module load gcc/6.2.0 boost/1.62.0 python/2.7.12 deeptools/3.0.2 && 
       bamCoverage --bam $input.bam  --binSize 25 
                   --normalizeUsing RPKM --outFileFormat bigwig
                   --scaleFactor 1 --numberOfProcessors 3
                   -o bigwigs/${basename}.bw
     ""","BamCoverage"
 }
 forward input
}

run {    
   branches * [Fastqc +
               Trimgalor +
               hisat2_map +
               sortAndIndex + 
               stringTie + 
               bamCoverage 

             ]
}



configfile: "config.yaml"

rule all: 
    input: 
        expand("{genome}.fasta", genome = config['GENOME']),
        expand("{genome}.fasta.fai", genome = config['GENOME']),
        expand("{chr}.vcf", chr =config['CHROMOSOMES']),
        expand("{COHORT}.hg38_multianno.vcf", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_filtered.vcf", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_sitesonly.vcf.gz", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_cohort_indels.recal", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_cohort_indels.tranches", COHORT=config['COHORT_NAME']),
        expand("{COHORT}.indel.recalibrated.vcf.gz", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_cohort_snps.recal", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_cohort_snps.tranches", COHORT=config['COHORT_NAME']),
        expand("{COHORT}.snps.recalibrated.vcf.gz", COHORT=config['COHORT_NAME'])


rule download: 
     params:
        config['GENOME']
     output:
        expand("{genome}.fasta", genome = config['GENOME']),
        config['DBSNP'],
        config['INDELS'],
        config['GOLD_STANDARD'],
        config['Axiom_Exome'],
        config['G_phase1'],
        config['G_omni2'],
        config['hapmap_3']

     shell:
         """
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
          wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
        """ 
rule index:
     input: 
           expand("{genome}.fasta", genome = config['GENOME']) 
     params: 
        config['GENOME'] 
     output:
        expand("{genome}.fasta.fai", genome = config['GENOME']),        
        expand("{genome}.fasta.rev.1.bt2", genome = config['GENOME']), 
        expand("{genome}.fasta.rev.2.bt2", genome = config['GENOME']),
        expand("{genome}.1.bt2", genome = config['GENOME']),
        expand("{genome}.2.bt2", genome = config['GENOME']),
        expand("{genome}.3.bt2", genome = config['GENOME']),
        expand("{genome}.4.bt2", genome = config['GENOME']), 

     shell: 
         """
          bowtie2-build {input} {params}
          samtools faidx {input} 
         """ 
rule trim: 
    input: 
       r1 = "{sample}.r_1.fq.gz",
       r2 = "{sample}.r_2.fq.gz"
    output: 
      val1 = "galore/{sample}.r_1_val_1.fq.gz",
      val2 = "galore/{sample}.r_2_val_2.fq.gz"
    shell: 
        """
         trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
        """ 

rule tosam:
    input:
        genome = expand("{genome}.fasta", genome=config['GENOME']),
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    output:
        '{sample}.sam'
    shell:
        "bowtie2 -x {input.genome} -1 {input.r1} -2 {input.r2} -S {output}"

rule AddRG: 
    input: 
       '{sample}.sam'
    output: 
       '{sample}.RG.sam' 
    params: 
        RG = config['RG']
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=@{params} RGSM={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={params}_{wildcards.sample} VALIDATION_STRINGENCY=SILENT" 


rule dedup: 
     input: 
         '{sample}.RG.sam'
     output:
       '{sample}.dedupped.bam',
       '{sample}.output.metrics'
     shell:
        "picard MarkDuplicates I={input} O={output[0]} CREATE_INDEX=true M={output[1]}"


rule recalibrate: 
    input:  
         sample = '{sample}.dedupped.bam', 
         genome = expand("{genome}.fasta", genome=config['GENOME'])
    output:
         '{sample}.recal_data.table',
         '{sample}.recalibrated.bam'
    params: 
        mem = "-Xmx800g",
        knownsites1 = config['DBSNP'],
        knownsites2 = config['INDELS'],
        knownsites3 = config['GOLD_STANDARD']
   
    shell:
       """
       gatk --java-options {params.mem} BaseRecalibrator -I {input.sample} -R {input.genome} --known-sites {params.knownsites1} --known-sites {params.knownsites2} --known-sites {params.knownsites3} -O {output[0]}
       gatk --java-options {params.mem} ApplyBQSR  -I {input.sample} -R {input.genome} --bqsr-recal-file {output[0]} -O {output[1]} 
       """ 

rule tovcf:
   input:
      expand('{sample}.recalibrated.bam',sample=config['SAMPLES']),
      genome = expand("{genome}.fasta", genome=config['GENOME'])

   params:
     mem_threads = {"-Xmx500g -XX:ParallelGCThreads=4"},
     I =  lambda w: "-I " + " -I ".join(expand("{sample}.recalibrated.bam", sample=config['SAMPLES'])),
   output:
       "{chr}.vcf" 
   shell:
       """
       gatk --java-options "{params.mem_threads}" HaplotypeCaller -R {input.genome} {params.I} -L {wildcards.chr}  -O {output[0]}
       """

rule merge:
     input: 
         expand("{chr}.vcf", chr = config['CHROMOSOMES'])  
     params: 
         mem = {"-Xmx600g"},
         I =  lambda w: "-I " + " -I ".join(expand("{chr}.vcf", chr =config['CHROMOSOMES']))
     output: 
         vcf = expand("{COHORT}.vcf", COHORT=config['COHORT_NAME'])
     shell: 
         """
         gatk --java-options "{params.mem}" MergeVcfs {params.I} -O {output}
         """ 
        
rule annotate: 
     input: 
         vcf = expand("{COHORT}.vcf", COHORT=config['COHORT_NAME'] )
     params: 
          humanDB = config['humanDB'],
          version = config['version'], 
          ANNOVAR = config['ANNOVAR'],
          output = config['COHORT_NAME'] 
     output:
         expand("{COHORT}.hg38_multianno.vcf", COHORT = config['COHORT_NAME'] )
     shell: 
          """
          {params.ANNOVAR}/table_annovar.pl {input.vcf} {params.humanDB} -buildver {params.version} -out {params.output} -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad_exome,avsnp147,dbnsfp33a,clinvar_20170130,revel -operation g,g,f,f,f,f,f,f,f  -nastring . -vcfinput
          """

rule hard_filter: 
    input: 
         vcf = expand("{COHORT}.hg38_multianno.vcf", COHORT = config['COHORT_NAME'])
    params:
        qd = config['QD'], 
        qual = config['QUAL'],
        sor = config['SOR'],
        fs = config['FS'], 
        mq = config['MQ'], 
        mqranksum = config['MQRankSum'],
        readposranksum = config['ReadPosRankSum']
    output:
         output = expand("{COHORT}_filtered.vcf", COHORT=config['COHORT_NAME'])
    shell: 
         """
         gatk VariantFiltration \
         -V {input} \
         -filter "QD < {params.qd}" --filter-name "QD2" \
         -filter "QUAL < {params.qual}" --filter-name "QUAL30" \
         -filter "SOR > {params.sor}" --filter-name "SOR3" \
         -filter "FS > {params.fs}" --filter-name "FS60" \
         -filter "MQ < {params.mq}" --filter-name "MQ40" \
         -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum-12.5" \
         -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum-8" \
         -O {output}
         """



rule hardfilter_ExcessHet: 
     input: 
          expand("{COHORT}_filtered.vcf", COHORT=config['COHORT_NAME'] )
     params: 
         ExcessHet = config['ExcessHet']
     output: 
          expand("{COHORT}_excesshet.vcf.gz", COHORT=config['COHORT_NAME'])             
     shell: 
        """
           gatk --java-options "-Xmx3g -Xms3g" VariantFiltration  -V {input} \
           --filter-expression "ExcessHet > {params}" \
           --filter-name ExcessHet -O {output} 
        """ 

rule sites_only: 
     input: 
        expand("{COHORT}_excesshet.vcf.gz", COHORT=config['COHORT_NAME'])
     output: 
        expand("{COHORT}_sitesonly.vcf.gz", COHORT=config['COHORT_NAME'])
     shell: 
       """ 
          gatk MakeSitesOnlyVcf -I {input} -O {output} 
       """

rule indels_tranches: 
    input: 
        expand("{COHORT}_sitesonly.vcf.gz", COHORT=config['COHORT_NAME'])
    output: 
        expand("{COHORT}_cohort_indels.recal", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_cohort_indels.tranches", COHORT=config['COHORT_NAME'])
    params: 
      DBSNP = config['DBSNP'], 
      GOLD_STANDARD = config['GOLD_STANDARD'],
      Axiom_Exome = config['Axiom_Exome'] 
    shell: 
      """   
       gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator -V {input} --trust-all-polymorphic \
       -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
       -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
       -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL --max-gaussians 4 \
       -resource:mills,known=false,training=true,truth=true,prior=12 {params.GOLD_STANDARD} \
       -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {params.Axiom_Exome} \
       -resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.DBSNP} -O {output[0]} --tranches-file {output[1]}
      """

rule SNPs_tranches:
    input:
        expand("{COHORT}_sitesonly.vcf.gz", COHORT=config['COHORT_NAME'])
    output: 
        expand("{COHORT}_cohort_snps.recal", COHORT=config['COHORT_NAME']),
        expand("{COHORT}_cohort_snps.tranches", COHORT=config['COHORT_NAME'])
    params:
        DBSNP = config['DBSNP'],
        G_phase1 = config['G_phase1'], 
        G_omni2 = config['G_omni2'],
        hapmap_3 = config['hapmap_3']  
    shell: 
       """ 
       gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator -V {input} --trust-all-polymorphic \
       -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
       -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
       -mode SNP --max-gaussians 6 \
       -resource:hapmap,known=false,training=true,truth=true,prior=15 {params.hapmap_3} \
       -resource:omni,known=false,training=true,truth=true,prior=12 {params.G_omni2} \
       -resource:1000G,known=false,training=true,truth=false,prior=10 {params.G_phase1} \
       -resource:dbsnp,known=true,training=false,truth=false,prior=7 {params.DBSNP} \
       -O {output[0]} --tranches-file {output[1]} 
      """


rule VQSLOD_indels: 
   input: 
      expand("{COHORT}_excesshet.vcf.gz", COHORT=config['COHORT_NAME']),
      expand("{COHORT}_cohort_indels.recal", COHORT=config['COHORT_NAME']),
      expand("{COHORT}_cohort_indels.tranches", COHORT=config['COHORT_NAME'])
   output: 
      expand("{COHORT}.indel.recalibrated.vcf.gz", COHORT=config['COHORT_NAME'])   
   shell:
    """ 
    gatk --java-options "-Xmx5g -Xms5g" \
    ApplyVQSR \
    -V {input[0]} \
    --recal-file {input[1]} \
    --tranches-file {input[2]} \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL \
    -O {output}
    """

rule VQSLOD_SNPs: 
   input: 
      expand("{COHORT}.indel.recalibrated.vcf.gz", COHORT=config['COHORT_NAME']),
      expand("{COHORT}_cohort_snps.recal", COHORT=config['COHORT_NAME']),
      expand("{COHORT}_cohort_snps.tranches", COHORT=config['COHORT_NAME'])
   output: 
      expand("{COHORT}.snps.recalibrated.vcf.gz", COHORT=config['COHORT_NAME'])
   shell: 
      """ 
      gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
      -V {input[0]} --recal-file {input[1]} --tranches-file {input[2]} \
      --truth-sensitivity-filter-level 99.7 \
      --create-output-variant-index true \
      -mode SNP \
      -O {output}  
      """ 

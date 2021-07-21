 #!/bin/bash

 
 
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -zxvf annovar.latest.tar.gz
cd annovar
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_exome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20170130 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar mitimpact2 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20170202 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar knownGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic70 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar nci60 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20170501 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp144 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp33a humandb/

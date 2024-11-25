REF = /nfs/turbo/umms-thahoang/sherine/zebrafish2/Drerio_genome 
FASTQ = /nfs/turbo/umms-thahoang/sherine/zebrafish2/FASTQs
S1=12222-ZF-1-GEX
S2=12222-ZF-2-GEX
CELLRANGER = /nfs/turbo/umms-thahoang/Tools/cellranger-9.0.0/bin/
Danio_rerio.GRCz11.105.gtf:
	wget http://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/Danio_rerio.GRCz11.105.gtf.gz
	gunzip Danio_rerio.GRCz11.105.gtf.gz


Danio_rerio.GRCz11.dna.primary_assembly.fa:
	wget http://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
	gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz


Danio_rerio.GRCz11.105.filtered.gtf:
	/nfs/turbo/umms-thahoang/sherine/tools/cellranger-7.2.0/cellranger mkgtf Danio_rerio.GRCz11.105.gtf Danio_rerio.GRCz11.105.filtered.gtf --attribute=gene_biotype:protein_coding


Drerio_genome:
	/nfs/turbo/umms-thahoang/sherine/tools/cellranger-7.2.0/cellranger mkref --genome=Drerio_genome --fasta=Danio_rerio.GRCz11.dna.primary_assembly.fa --genes=Danio_rerio.GRCz11.105.filtered.gtf

ZF1/outs:
	${CELLRANGER}/cellranger count   --id=ZF1 --transcriptome=${REF} --fastqs=${FASTQ}  --expect-cells=10000 --sample=${S1}  --create-bam  false 
	${CELLRANGER}/cellranger count   --id=ZF2 --transcriptome=${REF} --fastqs=${FASTQ}  --expect-cells=10000 --sample=${S2}  --create-bam  false    

ZF1_2:
	${CELLRANGER}/cellranger aggr --id=ZF1_2 --csv=aggr.csv --normalize=none




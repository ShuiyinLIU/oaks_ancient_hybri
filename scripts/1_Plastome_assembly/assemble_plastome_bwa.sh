# modified from https://github.com/ryanafolk/Assembly-tools/blob/master/BWA_assemble.sh

######################        # Assemble       #########################
for SLURM_ARRAY_TASK_ID in {1..4}; do

echo This is task $SLURM_ARRAY_TASK_ID

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

#Naming
R1FILE=$(ls reads/*_P1.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
name=`basename ${R1FILE} _P1.fastq`
echo "Running BWA on ${name}"
forward="reads/${name}_P1.fastq"
reverse="reads/${name}_P2.fastq"

### IR
# Assemble
# Make sure to index reference first
bwa mem -t 5 ../reference_LSR_SSR_IR/CYedithiae_KU382355.1_IR.fasta $forward $reverse > temp.sam

# Call consensus sequence
samtools sort -T ./tmp/aln.sorted -o temp.bam temp.sam
samtools mpileup -uf ../reference_LSR_SSR_IR/CYedithiae_KU382355.1_IR.fasta temp.bam > bcf.tmp
bcftools call --ploidy 1 -c bcf.tmp | vcfutils.pl vcf2fq > ./consensus_sequences/${name}.consensus.IR.fq
rm temp.sam temp.bam bcf.tmp ./tmp/aln.sorted*

### SSR
# Assemble
# Make sure to index reference first
bwa mem -t 5 ../reference_LSR_SSR_IR/CYedithiae_KU382355.1_SSR.fasta $forward $reverse > temp.sam

# Call consensus sequence
samtools sort -T ./tmp/aln.sorted -o temp.bam temp.sam
samtools mpileup -uf ../reference_LSR_SSR_IR/CYedithiae_KU382355.1_SSR.fasta temp.bam > bcf.tmp
bcftools call --ploidy 1 -c bcf.tmp | vcfutils.pl vcf2fq > ./consensus_sequences/${name}.consensus.SSR.fq
rm temp.sam temp.bam bcf.tmp ./tmp/aln.sorted*

### LSR
# Assemble
# Make sure to index reference first
bwa mem -t 5 ../reference_LSR_SSR_IR/CYedithiae_KU382355.1_LSR.fasta $forward $reverse > temp.sam

# Call consensus sequence
samtools sort -T ./tmp/aln.sorted -o temp.bam temp.sam
samtools mpileup -uf ../reference_LSR_SSR_IR/CYedithiae_KU382355.1_LSR.fasta temp.bam > bcf.tmp
bcftools call --ploidy 1 -c bcf.tmp | vcfutils.pl vcf2fq > ./consensus_sequences/${name}.consensus.LSR.fq
rm temp.sam temp.bam bcf.tmp ./tmp/aln.sorted*

done


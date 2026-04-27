module load Bowtie2

mkdir /gpfs/gibbs/project/noonan/yj345/makeHumanizedFasta/genome/hmm10
cd /gpfs/gibbs/project/noonan/yj345/makeHumanizedFasta/genome/hmm10

bowtie2-build hmm10.fa hmm10_index

mkdir /gpfs/gibbs/project/noonan/yj345/makeHumanizedFasta/genome/cmm10
cd /gpfs/gibbs/project/noonan/yj345/makeHumanizedFasta/genome/cmm10

bowtie2-build cmm10.fa cmm10_index

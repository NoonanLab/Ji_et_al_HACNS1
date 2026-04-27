conda activate sra

#Md
prefetch SRR4835044
fasterq-dump SRR4835044 -e 8 -p -O .

#PA2
prefetch SRR4835045
fasterq-dump SRR4835045 -e 8 -p -O .

#FL
prefetch SRR6294689
fasterq-dump SRR6294689 -e 8 -p -O .
prefetch SRR6294690
fasterq-dump SRR6294690 -e 8 -p -O .

#HL
prefetch SRR3950341
fasterq-dump SRR3950341 -e 8 -p -O .
prefetch SRR3950342
fasterq-dump SRR3950342 -e 8 -p -O .

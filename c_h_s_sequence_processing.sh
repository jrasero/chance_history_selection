#sequence data manipulation for chance, history, selection manuscript
nano allison_history.sh
#!/bin/bash
#SBATCH -p nodes
#SBATCH -J allison
#SBATCH -c 8

#run trimmomatic
module load trimmomatic/trimmomatic-0.36
for i in {01..48}; do trimmomatic PE -threads 8 -trimlog trim.log -phred33 ~/*_"$i"/*_R1_001.fastq.gz ~/_"$i"/*_R2_001.fastq.gz ~/trim/"$i"_forward_paired.fq.gz ~/trim/"$i"_forward_unpaired.fq.gz ~/trim/"$i"_reverse_paired.fq.gz ~/trim/"$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/cwm47/build/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done
#run breseq
module load breseq/breseq-0.31.0
for i in {01..48}; do time breseq -p -r ~/Abaumannii17978/NZ_CP012004_GCF_001077675.1_ASM107767v1_NC009084pAB1_009083pAB2_genomic.gbff  ~/trim/"$i"_forward_paired.fq.gz ~/trim/"$i"_forward_unpaired.fq.gz ~/trim/"$i"_reverse_paired.fq.gz -o ~/populations/breseq_"$i" -j 8; done 

#make sure full outputs for all samples
for i in {01..48}; do echo breseq_$i && ls breseq_$i/output/; done

~/gscripts/BreseqCat.py -p -d ~/populations

#subtract reference out
for i in {01..48}; do gdtools SUBTRACT -o ~/populations/"$i"_subtractref_output.gd ~/populations/breseq_"$i"/output/output.gd ~/abaum/uo1/alfonso_cipro_evol/breseq_atcc/NZ_CP012004_pAB123_ref/1_10p_v_vc17978/output/output.gd; done
#annotate each gdtools file
for i in {01..48}; do gdtools ANNOTATE -o "$i"_subtractref.txt -f TSV -r ~/Abaumannii17978/NZ_CP012004_GCF_001077675.1_ASM107767v1_NC009084pAB1_009083pAB2_genomic.gbff ~/populations/"$i"_subtractref_output.gd; done

#change outputs to sample name - run locally
for i in {01..48}; do sed -i '' "s/output/sample_$i/g" "$i"_subtractref.txt; done
#move analysis into R

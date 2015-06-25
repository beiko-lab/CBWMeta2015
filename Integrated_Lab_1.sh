#Canadian Bioinformatics Workshop 2015
#Marker Gene Integrated Assignment

#Make the script noisy
set -v

#First, we will link our sequence and reference data into our workspace
mkdir sequence_files
mkdir reference_data
mkdir scripts
ln -s ~/CourseData/integrated_assignment_day1/*.fastq.gz sequence_files/
ln -s ~/CourseData/integrated_assignment_day1/nwcs_arct_metadata.tsv .
ln -s ~/CourseData/integrated_assignment_day1/97* reference_data/
ln -s ~/CourseData/integrated_assignment_day1/core_set_aligned.fasta.imputed reference_data/

#Fetch some scripts from GitHub
wget -O scripts/mesas-pcoa https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-pcoa
wget -O scripts/mesas-uc2clust https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-uc2clust
#Mark the scripts as executable
chmod u+x scripts/*

#Assemble with PEAR
#Breakdown of the craziness below:
#  1) Find all gzipped FASTQ files and print their names
#  2) Use sed to only print out filename prefixes
#  3) Use sort to order the prefixes
#  4) Use uniq to remove duplicates (since fwd and rev files)
#  5) Take these files 4 at a time (since we have 4 cores)
#  6) For each file prefix, run PEAR in a subshell (multiprocessors)
#  7) Sleep for 60 seconds to let the 4 finish before starting next
#  8) Pipe the output to /dev/null so we don't get spammed
for i in "1" "5" "9" "13" "17" "21" "25"
do
    echo $i
    find ~/workspace/assignment1/sequence_files -name "*.fastq.gz" -printf '%f\n' | sed 's/_.*//' | sort | uniq | sed -n $i,$((i+3))p | while read line; do ( pear -f sequence_files/${line}_1.fastq.gz -r sequence_files/${line}_2.fastq.gz -o ${line} & ); done > /dev/null
    sleep 60
done

#Give a wee bit more time just in case any of the PEAR jobs have not caught up
#Sleep until we have
while [ $( ps -ef | grep pear | wc -l ) -gt 1 ]
do
 sleep 10
 echo "Waiting on PEAR to finish..."
done

#Ensure the sequence file is no longer here
rm -rf seq.fasta

#Take in the FASTQ files, funnel them into a single FASTA file
#FASTA headers in the form: ">SAMPLEID_SEQUENCENUMBER"
for filename in $( ls *.assembled.fastq )
do
    awk 'BEGIN{ORS=""; i=0;}{split(FILENAME, x, "."); prefix=x[1]; sub("@","",prefix); print ">" prefix "_" i "\n"; i+=1; getline; print; print "\n"; getline; getline;}' ${filename} >> seq.fasta
done

#Clean up our directory
mv *.fastq sequence_files

#Cluster with UPARSE pipeline via USEARCH

#Collapse dataset to remove redundancy
usearch -derep_fulllength seq.fasta -fastaout derep.fa -sizeout
#Sort the sequences by abundance
usearch -sortbysize derep.fa -fastaout sorted.fa -minsize 2
#Cluster the sorted, dereplicated, minsize=2 sequence set
usearch -cluster_otus sorted.fa -otus otus.fa -otu_radius_pct 3 -sizeout -uparseout results.txt
#Rename the OTUs for QIIME
awk 'BEGIN{count=0;}{if ($0~/>/){print ">" count; count+=1;} else {print}}' otus.fa > rep_set.fasta
#Map all sequences (even singletons) back onto the OTUs
usearch -usearch_global seq.fasta -db rep_set.fasta -strand both -id 0.97 -uc map.uc -threads 4
#Convert USEARCH output to QIIME OTU list format
scripts/mesas-uc2clust -t 4 map.uc seq_otus.txt

#Clean up our directory
mkdir cluster
mv results.txt otus.fa map.uc rep_set.fasta seq_otus.txt sorted.fa derep.fa cluster

#Switch to use the right BIOM package
export PYTHONPATH=~/local/lib/python2.7/site-packages
#Tell QIIME where RDP is
export RDP_JAR_PATH=/usr/local/rdp_classifier_2.2/rdp_classifier-2.2.jar

#QIIME commands

#Assigning taxonomic classifications via RDP v2.2
assign_taxonomy.py -m rdp -i cluster/rep_set.fasta -o taxonomic_classification -t reference_data/97_otu_taxonomy.txt -r reference_data/97_otus.fasta -c 0.6
#Align sequences to a template for tree generation
align_seqs.py -m pynast -i cluster/rep_set.fasta -o alignment -t reference_data/core_set_aligned.fasta.imputed
#Filter empty columns out of alignment
filter_alignment.py -i alignment/rep_set_aligned.fasta -o alignment -s
#Make a directory for the tree
mkdir phylogeny
make_phylogeny.py -i alignment/rep_set_aligned_pfiltered.fasta -t fasttree -o phylogeny/rep_set.tre -l phylogeny/log.txt
#Make OTU table containing counts and taxonomic classifications
mkdir otu_table
make_otu_table.py -i cluster/seq_otus.txt -o otu_table/otu_table.biom -t taxonomic_classification/rep_set_tax_assignments.txt
#Convert the BIOM format file to tab-separated format (human readable)
biom convert --to-tsv -i otu_table/otu_table.biom --table-type='OTU table' -o otu_table/otu_table.tab --header-key=taxonomy --output-metadata-id=Consensus\ Lineage
#We use awk on the converted OTU table to determine the lowest sequence depth
#This is passed as a parameter to QIIME's rarefaction script
single_rarefaction.py -i otu_table/otu_table.biom -o otu_table/otu_table_rarefied.biom -d `awk 'BEGIN{FS="\t"} NR == 1 { } NR == 2 { max = NF-1; } NR > 2 { for (i = 2; i <= max; i++) { c[i] += $i; } } END { smallest = c[2]; for (i = 3; i <= max; i++) { if (c[i] < smallest) { smallest = c[i]; }} print smallest; }' otu_table/otu_table.tab`
#Convert rarefied table to tab-separated format
biom convert --to-tsv -i otu_table/otu_table_rarefied.biom --table-type='OTU table' -o otu_table/otu_table_rarefied.tab --header-key=taxonomy --output-metadata-id=Consensus\ Lineage

#We now have OTU tables! We can do analyses on it now.

#Taxonomy bar plots via QIIME
#summarize_taxa_through_plots.py -f -s -i otu_table/otu_table.biom -o taxaplot -m nwcs_arct_metadata.tsv

#Create distance matrices
beta_diversity.py -i otu_table/otu_table_rarefied.biom -o qiime_pcoa/distance/ -m weighted_unifrac -t phylogeny/rep_set.tre
beta_diversity.py -i otu_table/otu_table_rarefied.biom -o qiime_pcoa/distance/ -m bray_curtis

#Run PCoA on these matrices
principal_coordinates.py -i qiime_pcoa/distance/ -o qiime_pcoa/pcoa/

#Plot them
make_2d_plots.py -i qiime_pcoa/pcoa -m nwcs_arct_metadata.tsv -o qiime_pcoa/pcoa_plots/
chmod -R o+x qiime_pcoa/pcoa_plots

#Run an R script for PCoA
mkdir mesas_pcoa
scripts/mesas-pcoa -i otu_table/otu_table_rarefied.tab -m nwcs_arct_metadata.tsv -d bray -o mesas_pcoa/pcoa-bray.pdf


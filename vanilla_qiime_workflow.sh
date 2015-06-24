#"Vanilla" QIIME workflow

#What you need to have done before starting:
#Demultiplexed, assembled, quality checked, and filtered your sequences
#Gathered all your metadata into a tab-separated spreadsheet file (.tsv)

#Required - FASTA file with assembled/quality filtered sequences, 
#Files:      and FASTA headers with format ">SampleID_SequenceNumber"
#            (ie, ">SampleA_143")
#          - Metadata mapping file (tab separated) with first line as 
#            header with "#SampleID" as first column, containing sample
#            IDs corresponding to FASTA file

#All scripts documented at http://qiime.org/scripts

#Let's call our sequence file "seq.fasta" and mapping file "mapping.tsv"

mkdir cluster
#Use the newest USEARCH QIIME offers, 97% identity (-s 0.97), 
#discard singletons for clustering (-g 2), check the reverse complement (-z),
#and select seeds starting with most abundant
#Please note that this will not work in our cloud session since usearch 6.1 is not installed
pick_otus.py -i seq.fasta -m usearch61 -o cluster -s 0.97 -g 2 -z --sizeorder

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

#At this point you need to check what the lowest sampling depth is
#This can be determined by summarizing the OTU table and inspecting the
#output file (otu_table/otu_table_summary.txt)
biom summarize-table -i otu_table/otu_table.biom -o otu_table/otu_table_summary.txt

#After checking the summary table for the minimum number of sequences in
#a sample, replace the value in the command below
single_rarefaction.py -i otu_table/otu_table.biom -o otu_table/otu_table_rarefied.biom -d <INSERT DEPTH HERE>

#Create taxonomic bar plots
summarize_taxa_through_plots.py -f -s -i otu_table/otu_table.biom -o taxaplot -m mapping.tsv

#Create distance matrix for weighted UniFrac
beta_diversity.py -i otu_table/otu_table_rarefied.biom -o qiime_pcoa/distance/ -m weighted_unifrac -t phylogeny/rep_set.tre

#Run PCoA on these matrices
principal_coordinates.py -i qiime_pcoa/distance/ -o qiime_pcoa/pcoa/

#Plot them
make_2d_plots.py -i qiime_pcoa/pcoa -m mapping.tsv -o qiime_pcoa/pcoa_plots/

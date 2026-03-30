# Part_II_Project
#Submitted on 30/3/2026
#Workflow for curating a universal TE library and extracting high-confidence consensuses

# Notes:
#All commands are for Onix/macOS unless otherwise stated
#Ensure all required scripts (libstats.R, Filter_libstats.py, passing_clusters.py, split_clusters_by_source.py, full_TE_pipeline.sh, Generate_consensus.R) are in the working directory
#Tools: CD-HIT, R, Python3, TE-Aid, EMBOSS, RepeatMasker, MAFFT


# Step 1: Curate a universal library from the 4 source libraries


# Source libraries
#FlyBase: Berg.fa
#FlyBase D. mel: Berg_Dmel.fa
#MCH: MCH.fa
#Pantera: Pantera.fa
#Gonzalez: Pepi.fa

# 1.1 Combine Berg and Berg_Dmel to assess redundancy
cat Berg.fa Berg_Dmel.fa > FlyBase_combined.fa

# Cluster FlyBase combined library using 80-80 rule
cd-hit-est -i FlyBase_combined.fa -o FlyBase_Clustered_80_80 -c 0.8 -aS 0.8 -d 0
#Output: FlyBase_Clustered_80_80.fa and FlyBase_Clustered_80_80.fa.clstr

# 1.2 Add prefixes to each library (to keep track of source)
awk '/^>/ {print ">" "berg_" substr($0,2); next} {print}' Berg.fa > Berg_prefixed.fa
awk '/^>/ {print ">" "mch_" substr($0,2); next} {print}' MCH.fa > MCH_prefixed.fa
awk '/^>/ {print ">" "pantera_" substr($0,2); next} {print}' Pantera.fa > Pantera_prefixed.fa
awk '/^>/ {print ">" "pepi_" substr($0,2); next} {print}' Pepi.fa > Pepi_prefixed.fa


# Step 2: Cluster the universal library

cat Berg_prefixed.fa MCH_prefixed.fa Pantera_prefixed.fa Pepi_prefixed.fa > universal_prefixed.fa

cd-hit-est -i universal_prefixed.fa -o Universal_prefixed_Clustered_80_80 -c 0.8 -aS 0.8 -d 0
#Output: Universal_prefixed_Clustered_80_80.fa and .clstr


# Step 3: Filter for high-confidence TE entries

# Run libstats.R for each library to generate stats
Rscript libstats.R Berg_prefixed.fa
Rscript libstats.R MCH_prefixed.fa
Rscript libstats.R Pantera_prefixed.fa
Rscript libstats.R Pepi_prefixed.fa

# Filter libraries using Python script (Filter_libstats.py)
python3 Filter_libstats.py Berg_prefixed.fa.stats
python3 Filter_libstats.py MCH_prefixed.fa.stats
python3 Filter_libstats.py Pantera_prefixed.fa.stats
python3 Filter_libstats.py Pepi_prefixed.fa.stats
#Output example: Berg_prefixed.fa.stats.pass / Berg_prefixed.fa.stats.fail


# Step 4: Map passed/failed entries back to clusters

cat Berg_prefixed.fa.stats.pass MCH_prefixed.fa.stats.pass Pantera_prefixed.fa.stats.pass Pepi_prefixed.fa.stats.pass > all_pass.stats
cat Berg_prefixed.fa.stats.fail MCH_prefixed.fa.stats.fail Pantera_prefixed.fa.stats.fail Pepi_prefixed.fa.stats.fail > all_fail.stats

# Map to clusters
python3 passing_clusters.py all_pass.stats all_fail.stats Universal_prefixed_Clustered_80_80.fa.clstr
# Output: clusters_with_pass.tsv and clusters_without_pass.tsv


# Step 4.2: Select one consensus per family

# Semi-automatic prioritization based on source and quality
python3 split_clusters_by_source.py
#Outputs:
#1) clusters_with_pass_mch_pantera_only.tsv 
#2) clusters_with_pass_pepi_berg_only.tsv
#3) clusters_with_pass_mixed.tsv

#Post-processing:
#- Discard fail / low-quality entries as described
#- Keep only entries of interest

# Merge all retained clusters
cat clusters_with_pass_mch_pantera_only.tsv clusters_with_pass_pepi_berg_only.tsv clusters_with_pass_mixed.tsv > universal_passed_entries.tsv


# Step 4.3: Extract names of passed entries

cut -f2 universal_passed_entries.tsv > passed_names.txt


# Step 4.4: Filter universal FASTA library using conserved FASTA names
awk 'BEGIN{while(getline<"passed_names.txt") names[$1]=1} 
/^>/{f=names[substr($0,2)]} f{print}' universal.fa > universal_passed.fa
#Output: universal_passed.fa (final high-confidence universal TE library)

# Step 5: Masking against the genome
RepeatMasker -pa 8 -lib universal_passed.fa -dir RM_output Dame_genome.fa
Output:RM_output/Dame_genome.fa.out
./full_TE_pipeline.sh RM_output/Dame_genome.fa.out universal_passed.fa Dame_genome.fa D.americana Dame

# Step 6: Manual validation with TE-Aid
TE-aid analyze -i Dame_Pipeline/consensus/consensi.fa  -o TE-aid_results -s D_americana

# Step 7: Phylogenetic Tree
grep -A 1 "LINE" universal_library.fa > LINE_elements.fasta
getorf -sequence LINE_elements.fasta -outseq LINE_ORFs.fasta -minsize 2000
mafft --auto LINE_ORFs.fasta > LINE_ORF2_aln.fasta
iqtree -s LINE_ORF2_aln.fasta -st AA -m MFP -bb 1000 -alrt 1000 -nt AUTO

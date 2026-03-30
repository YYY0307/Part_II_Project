# Part_II_Project
# Submitted on 30/3/2026

#Read relevant script/command according to Methodologies
#All comment in onix (macOS unless otherwise stated)
#Step 1: Curating a universal library from the 4 libraries

#Library naming: 
FlyBase Library: Berg.fa
FlyBase(D_mel Library): Berg_Dmel.fa
MCH Library: MCH.fa
pantera Library: pantera.fa
Gonzalez library: Pepi.fa

#1.1 Confirming the robustness of Berg_Dmel.fa against Berg.fa: 
cat Berg.fa Berg_dmel.fa > FlyBase_combined.fa 
cd-hit-est -i FlyBase_combined.fa -o FlyBase_Clustered_80_80 -c 0.8 -aS 0.8 -d 0 # According to the 80-80 rule, retention of the fasta header, as shown in the report

#Step 2: Filtering for high-confidence TE Entries
python3 libstats.R (library.fa)# name of the library (e.g. Berg.fa) #libstats.R attached in here


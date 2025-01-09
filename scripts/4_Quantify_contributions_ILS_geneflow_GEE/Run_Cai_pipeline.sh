#This script is used for Quantifying the relative importance of gene tree estimation error, ILS, and gene flow
#primarily following https://github.com/lmcai/Coalescent_simulation_and_gene_flow_detection

##############################################################################
### for dataset HYB-98RT
mkdir CS_FGdetection && mkdir CS_FGdetection/RT431_98 && cd CS_FGdetection/RT431_98


####  1\The dependent variable — Gene tree variation
#It quantifies gene tree topological variation across the tree.
while read name;do cat ${name} >>BSgenetree.trees;done <../../ASTRAL_bootstrap/RT431_98.bs-files
bsub -J TQ_IQgcf -q Q104C512G_X4  -R "span[hosts=1]" -o output.%J -e error.%J -n 1 iqtree -t astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr --gcf BSgenetree.trees --prefix concord
#output:
#(1)concord.cf.tree: Newick tree with gCF assigned for each internal branch of the reference tree. If the reference tree already has some 
#branch label (such as bootstrap support in this case), gCF will be appended to the existing label separated by a /.
#(2)concord.cf.branch: Newick tree with internal branch IDs.
#(3)concord.cf.stat: A tab-separated table with gCF and gDF (gene discordance factor) for every internal branch (rows of the table). The 
#ID column can be linked with concord.cf.branch file. This file can be read in R to do some plot (see below).



####  2\Gene tree estimation error
#It quantifies the anticipated level of analytical error for each node in the species tree, which is largely affected by the length of the 
#internal branch.
#renaming RAxML_info.*, RAxML_bipartitionsBranchLabels.*, RAxML_bestTree.*
while read name;do mv RAxML_info.${name}.raxml_bs.tre ${name}.info;
mv RAxML_bipartitionsBranchLabels.${name}.raxml_bs.tre ${name}.bipaBL.tre;
mv RAxML_bestTree.${name}.raxml_bs.tre ${name}.raxml_best.tre;done <../genelist_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.txt

#use Seq-Gen to simulate gene sequences based on species tree
#the input tree of Seq-Gen can not include node labels
cp astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr.noNlabel
sed -i "s/)[^:]*:/):/g" astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr.noNlabel
mkdir Seq_Gen

#simulate gene sequences
#RT431-98 (median:985.5, we use 986)
#Time: <1min
while read name;do echo ${name}
alpha_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "alpha: "|cut -f 2 -d " "`
alpha=`echo -n $alpha_raw`
rAC_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "rate A <-> C: "|cut -f 5 -d " "`
rAC=`echo -n $rAC_raw`
rAG_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "rate A <-> G: "|cut -f 5 -d " "`
rAG=`echo -n $rAG_raw`
rAT_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "rate A <-> T: "|cut -f 5 -d " "`
rAT=`echo -n $rAT_raw`
rCG_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "rate C <-> G: "|cut -f 5 -d " "`
rCG=`echo -n $rCG_raw`
rCT_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "rate C <-> T: "|cut -f 5 -d " "`
rCT=`echo -n $rCT_raw`
rGT_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "rate G <-> T: "|cut -f 5 -d " "`
rGT=`echo -n $rGT_raw`
fA_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "freq pi(A): "|cut -f 3 -d " "`
fA=`echo -n $fA_raw`
fC_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "freq pi(C): "|cut -f 3 -d " "`
fC=`echo -n $fC_raw`
fG_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "freq pi(G): "|cut -f 3 -d " "`
fG=`echo -n $fG_raw`
fT_raw=`cat ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.info|grep "freq pi(T): "|cut -f 3 -d " "`
fT=`echo -n $fT_raw`
~/applications/Seq-Gen-1.3.4/source/seq-gen -l 986 -m GTR -a ${alpha} -r ${rAC} ${rAG} ${rAT} ${rCG} ${rCT} ${rGT} -f ${fA} ${fC} ${fG} ${fT} -or < astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr.noNlabel >Seq_Gen/${name}.phy
done <../../genelist_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.txt

#run gene tree reconstrunction with raxml
cd Seq_Gen
bsub -J LSY_err -q Q104C512G_X4 -R "span[hosts=1]" -o output.%J -e error.%J -n 50 parallel -j 50 --eta raxmlHPC-PTHREADS-AVX -s {}.phy -n {}.raxml_bs.tre -m GTRGAMMA -f a -x 12345 -p 12345 -# 100 -T 2 :::: ../genelist_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.txt
cd ../

#calculate how often each node in the species tree is recovered in the inferred gene trees sim_gene.trees using the bipartition in RAxML.
cat Seq_Gen/RAxML_bipartitions.*>sim_gene.trees
bsub -J LSY_err -q Q104C512G_X4 -R "span[hosts=1]" -o output.%J -e error.%J -n 1 raxmlHPC -f b -t astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -z sim_gene.trees -m GTRGAMMA -n ERR -o JugCarya_pallida_P026_WE07,JugJuglans_australis_P026_WE05,MyrMorella_pringlei_P022_WH09,MyrMyrica_hartwegii_P022_WG03,CSmonbeigiana_trans,BetBetula_litwinowii_P021_WD06,BetBetula_pamirica_P021_WD05,NotNothofagus_fusca_P022_WG06



####  3\ILS
#It quantifies the level of incomplete lineage sorting in each node of the species tree. 
#remove all branch length and node label in astral tree, and then with this topology as the input of RAxML
bsub -J LSY -q Q104C512G_X4  -R "span[hosts=1]" -o output.%J -e error.%J -n 1 Rscript ~/liushuiyin/applications/Scripts_lsy/1.converse_tree_topology.R astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr topology_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr
cp /ds3200_1/users_root/yitingshuang/liushuiyin/phylogeneticsignal/tre/conRT_ot8ot7ot2R2.V3sra.tm_431.V1_98.phy ./
bsub -J LSY98 -q Q104C512G_X4  -R "span[hosts=1]" -o output.%J -e error.%J -n 5 raxmlHPC-PTHREADS-AVX2 -g topology_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -m GTRGAMMA -p 12345 -s conRT_ot8ot7ot2R2.V3sra.tm_431.V1_98.phy -n astralconstrained.tre -T 5
pxrr -t RAxML_bestTree.astralconstrained.tre -o RAxML_bestTree.astralconstrained.tre.rr -g JugCarya_pallida_P026_WE07,JugJuglans_australis_P026_WE05,MyrMorella_pringlei_P022_WH09,MyrMyrica_hartwegii_P022_WG03,CSmonbeigiana_trans,BetBetula_litwinowii_P021_WD06,BetBetula_pamirica_P021_WD05,NotNothofagus_fusca_P022_WG06



####  4\Gene flow
mkdir Geneflow && cd Geneflow

#(1) Simulate gene trees under coalescent model for each BP species tree
#The R package phybase we used for simulation do have quite stringent format requirement: rooted, no node labels, and 
#branch length has to be larger than 0.
printf "Simulating gene trees under MSC model...\n"
mkdir geneTr_sim
#assign branch length for terminal tips
#"add-bl.py" is from https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py
python ~/applications/Scripts_lsy/add-bl.py - astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre > astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.tre
python ~/applications/Scripts_lsy/add-bl.py - RT431_98.bootstrapped.100.astral.tre > RT431_98.bootstrapped.100.astral.add-bl.tre
#drop node labels
cp astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.tre astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre
sed -i "s/)[^:]*:/):/g" astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre
cp RT431_98.bootstrapped.100.astral.add-bl.tre RT431_98.bootstrapped.100.astral.add-bl.noNlabel.tre
sed -i "s/)[^:]*:/):/g" RT431_98.bootstrapped.100.astral.add-bl.noNlabel.tre
#reroot astral tree
pxrr -t astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre -o astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr -g NotNothofagus_fusca_P022_WG06 
pxrr -t RT431_98.bootstrapped.100.astral.add-bl.noNlabel.tre -o RT431_98.bootstrapped.100.astral.add-bl.noNlabel.tre.rr -g NotNothofagus_fusca_P022_WG06 
head -n 100 RT431_98.bootstrapped.100.astral.add-bl.noNlabel.tre.rr > 4BP.sp.rooted.trees
#combine all rooted empirical gene trees into a single file
cat ../Hyb_DB/genetrees_rr_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98/* > RT431_98/4ML.rooted.gene.trees
#run MSC_geneTr_simulator.R
Rscript --vanilla ~/applications/Coalescent_simulation_and_gene_flow_detection/MSC_geneTr_simulator.R astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr 4BP.sp.rooted.trees 4ML.rooted.gene.trees
rm geneTr_sim/*.tem.genetrees

#(2) Count triplet frequency in empirical gene trees and simulated gene trees
printf "Counting triplet frequency in the empirical data...\n"
python ~/applications/Coalescent_simulation_and_gene_flow_detection/triple_frequency_counter.py 4ML.rooted.gene.trees astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr
#excepted result: *.trp.tsv
#format: column1--species names of the triplet, sorted alphabetically; column2--triplet frequencies of (sp1,sp2);column3--triplet frequencies of (sp1,sp3);column4--triplet frequencies of (sp2,sp3).
printf "Counting triplet frequency in the simulated gene trees...\n"
ls geneTr_sim>list_geneTr_sim.txt
bsub -J LSY_tri89 -q Q104C512G_X4 -R "span[hosts=1]" -o output.%J -e error.%J -n 50 parallel -j 50 --eta python ~/applications/Coalescent_simulation_and_gene_flow_detection/triple_frequency_counter.py geneTr_sim/{} astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr :::: list_geneTr_sim.txt

#(3) Find significantly unbalanced triplets and map to species tree
mv astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr 4ML.astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr
bsub -J LSY_unbal -q Q104C512G_X4 -R "span[hosts=1]" -o output.%J -e error.%J -n 1 python ~/applications/Coalescent_simulation_and_gene_flow_detection/find_unbalanced_triplets.py 4ML.astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr
bsub -J LSY_map -q Q104C512G_X4 -R "span[hosts=1]" -o output.%J -e error.%J -n 1 python ~/applications/Coalescent_simulation_and_gene_flow_detection/triplet_mapper.py 4ML.astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.add-bl.noNlabel.tre.rr unbalanced.trp.tsv


######        Regression analysis in R         #######
#We modified the original R script "ralaimpo.R" to format the data above into a matrix (Make sure each row corresponds to the same internal node),
#and conduct importance decomposition.
#see "2.regression_analysis_relaimpo.R"



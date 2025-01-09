# This script is used to conduct coalescent simulations for characterizing ILS expectations

## rescaling with a scale factor of 4
pxtscale -t pruned215rl_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -s 4 -o pruned215rl.4x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr

## rescaling with a scale factor of 2
pxtscale -t pruned215rl_astral_MO_ot8ot7ot2R2.V3sra.tm_431.V1_89.tre.rr -s 2 -o pruned215rl.2x_astral_MO_ot8ot7ot2R2.V3sra.tm_431.V1_89.tre.rr

## run coalescent simulations in Dendropy
#"generateCoalescentTrees_preservetaxonnames.py" is from https://github.com/ryanafolk/tree_utilities
python ~/applications/simulate_gene_trees/generateCoalescentTrees_preservetaxonnames.py pruned215rl.4x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr 1000 genetrees_pruned215rl.4x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr
python ~/applications/simulate_gene_trees/generateCoalescentTrees_preservetaxonnames.py pruned215rl.2x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr 1000 genetrees_pruned215rl.2x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr


## replace all "[&R] " and "'" in two tree files with 1000 simulated gene trees in NotePad+.

## map simulated gene trees to nuclear and chloroplast species trees with RAxML
# on plastome tree
raxmlHPC -f b -t pruned215_RAxML_bipartitions.mp80_overlap324.ActionV2_4.223.sra_p0.3part.tre -z genetrees_pruned215rl.4x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -m GTRGAMMA -n mapCS.astralRT215_98.x4_RAxML_bipartitions.mp80_overlap324.ActionV2_4.223.sra_p0.3part.tre -o Carpinus_monbeigiana_NC_039997.1
raxmlHPC -f b -t pruned215_RAxML_bipartitions.mp80_overlap324.ActionV2_4.223.sra_p0.3part.tre -z genetrees_pruned215rl.2x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -m GTRGAMMA -n mapCS.astralRT215_98.x2_RAxML_bipartitions.mp80_overlap324.ActionV2_4.223.sra_p0.3part.tre -o Carpinus_monbeigiana_NC_039997.1

# on astral tree
raxmlHPC -f b -t pruned215rl_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -z genetrees_pruned215rl.4x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -m GTRGAMMA -n mapCS.astralRT215_98.x4_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -o Carpinus_monbeigiana_NC_039997.1
raxmlHPC -f b -t pruned215rl_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -z genetrees_pruned215rl.2x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -m GTRGAMMA -n mapCS.astralRT215_98.x2_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -o Carpinus_monbeigiana_NC_039997.1

#on con-ml tree
raxmlHPC -f b -t pruned215rl_RAxML_bipartitions.conRT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre -z genetrees_pruned215rl.4x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -m GTRGAMMA -n mapCS.astralRT215_98.x4_RAxML_bipartitions.conRT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre -o Carpinus_monbeigiana_NC_039997.1
raxmlHPC -f b -t pruned215rl_RAxML_bipartitions.conRT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre -z genetrees_pruned215rl.2x_astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr -m GTRGAMMA -n mapCS.astralRT215_98.x2_RAxML_bipartitions.conRT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre -o Carpinus_monbeigiana_NC_039997.1

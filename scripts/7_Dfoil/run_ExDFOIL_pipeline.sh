#################################################################################################################
#########         Dfoil analysis
## This script is primarily following https://github.com/SheaML/ExDFOIL


### for HYB-98RT
#####(1) Stage 1: Selection of Taxa
Rscript /data/liushuiyin/applications/ExDFOIL/Scripts/DFOIL_Picker.R namelist_dfoil_HYB_B.txt mcctree.calib0a9f_conMO431_89k_top20_10000resample.newick
mv myoutput.txt DFOIL_picked_B.txt #5 323 154 (HYB-RT98-B)

#filter 4-taxa combinations, which belong to the same group (i.e.,genus or section)
#5 316 699 (HYB-RT98-B); #224 998 for RNA-DB
python /data/liushuiyin/applications/Scripts_lsy/1.dfoil_filter_combinations.py DFOIL_picked_B.txt sampleinfo_HYB_B.txt DFOIL_picked_B_filtered.txt



#####(2) Stage 2: Running DFOIL with GNU Parallel
# 1. Getting site-pattern counts
#"2.fasta2foiler.sh" is modified from "https://github.com/SheaML/ExDFOIL/blob/master/Scripts/fasta2foiler.sh" by adding an additional paramater "cpu"
bsub -J TQ_dfoilaa -q Q104C512G_X4  -R "span[hosts=1]" -o output.%J -e error.%J -n 60 2.fasta2foiler.sh DFOIL_picked_B_filtered.txt conRT_ot8ot7ot2R2.V3sra.tm_431.V1_98.fa TR_OG.txt 60
#OUTPUTS: saved in folder "counts", each file in "counts" is formated as:
#chrom	      position	AAAAA	AAABA	AABAA	AABBA	ABAAA	ABABA	ABBAA	ABBBA	BAAAA	BAABA	BABAA	BABBA	BBAAA	BBABA	BBBAA	BBBBA
#/dev/fd/63	      0     27969    302	377	      77	365	     21	      13	  25	 446	  18	  14	  16	  203	  54	  20	 626

# 2. Running DFOIL
#"3.dfoiler_alt.sh" is modified from "https://github.com/SheaML/ExDFOIL/blob/master/Scripts/dfoiler_alt.sh" by adding an additional paramater "cpu"
bsub -J TQ_dfoil -q Q104C512G_X4  -R "span[hosts=1]" -o output.%J -e error.%J -n 35 3.dfoiler_alt.sh 35
#OUTPUTS：in the directory ./dfoil/ with the extension .dfoil_alt
#The "precheck" files are written to the folder ./precheck/ with the extension .precheck_alt
#e.g., oberon10_JJW683p1.oberonS1_SML142.ornatus3_JJW673.ornatusS1_SML128.dfoil_alt：
#chrom	coord	total	dtotal	T12	T34	T1234	DFO_left	DFO_right	DFO_total	DFO_stat	DFO_chisq	DFO_Pvalue	DIL_left	DIL_right	DIL_total	DIL_stat	DIL_chisq	DIL_Pvalue	DFI_left	DFI_right	DFI_total	DFI_stat	DFI_chisq	DFI_Pvalue	DOL_left	DOL_right	DOL_total	DOL_stat	DOL_chisq	DOL_Pvalue	introgression	introgna	intrognone	introg13	introg14	introg23	introg24	introg31	introg41	introg32	introg42	introg123	introg124
#/dev/fd/63	0	30494	97	0.007116154	0.014363481	0.017593625	30	32	62	-0.032258065	0.064516129	0.799495362	34	28	62	0.096774194	0.580645161	0.446059549	25	40	65	-0.230769231	3.461538462	0.062811848	29	36	65	-0.107692308	0.753846154	0.385261241	none	0	1	0	0	0	0	0	0	0	0	0	0

# 3. Summarizing tests
ls counts/ > counts_file.txt
awk '{print "counts/" $0}' counts_file.txt >counts_absfile.txt
bsub -J TQ_dfoil -q Q104C512G_X4  -o output.%J -e error.%J -n 45 parallel -j 45 'echo -n {}|summarizer_alt.sh' :::: counts_absfile.txt
#OUTPUTS: in a directory summary/ with the extension .summary_alt.
#To to get all the test summaries in a single file for use in the next stage:
find summary/ -name "*.summary_alt"|xargs -i cat {}>> all.summary.alt.txt
wc -l all.summary.alt.txt
head -n 2 all.summary.alt.txt



####(3) Stage 3: Summarizing and Visualizing Results
# 1. Associating Sample Info
bsub -J TQ_dfoil -q Q104C512G_X4  -o output.%J -e error.%J -n 1 associate.sh all.summary.alt.txt sampleinfo_HYB_B_space.txt
#OUTPUTS：all.summary_appended.txt
# 2. Getting DFOIL results into R and visualize
#see "4.visualize_DFOILres.R"




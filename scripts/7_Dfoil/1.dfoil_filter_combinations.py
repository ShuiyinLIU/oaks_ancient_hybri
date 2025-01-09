
#
import sys

def filter_combinations_withinGroup(dfoil_comb_f=sys.argv[1],sampleinfo_f=sys.argv[2],output_f=sys.argv[3]):
    #(1)creat a dict for sample whcih belongs to a certain group
    sample_dict={}
    with open(sampleinfo_f,'r') as sample_f:
        for sp in sample_f.readlines():
            taxa=sp.rstrip().split('\t')[0]
            #the third column is the sectionID info.
            sample_dict[taxa]=sp.rstrip().split('\t')[2]
    with open(output_f,'a+') as outf:
        with open(dfoil_comb_f,'r') as combf:
            n=0
            for comb in combf.readlines():
                n+=1
                print(n)
                group=[]
                #creat a list for 4 taxa on their group info.
                group.append(sample_dict[comb.rstrip().split(' ')[0]])
                group.append(sample_dict[comb.rstrip().split(' ')[1]])
                group.append(sample_dict[comb.rstrip().split(' ')[2]])
                group.append(sample_dict[comb.rstrip().split(' ')[3]])
                if len(set(group))>1:
                    outf.write(comb)
                else:
                    print(comb)

filter_combinations_withinGroup()








#written at 17:29, TUE, 2019-03-07 with Python3.6
# csvfile, can't include chinese words!!

#union
#the parameter 'csvfile' should contain path.
import sys
import os

def split_csvFile(csvfile=sys.argv[1], outputPath=sys.argv[2]):
    with open(csvfile,'r') as inputfile:
        sp_list = []
        #header=""
        for line in inputfile.readlines():
            if line.rstrip().split(',')[0] == 'species':
                #header = line
                pass
            else:
                sp = line.rstrip().split(',')[0]
                if sp not in sp_list:
                    sp_list.append(sp)
                else:
                    pass
        print("Toatally " + str(len(sp_list)) + " species!")
        #print(header)
        # a loop for writing all records of same species to same csv file
        for wuz in sp_list:
            print(wuz)
            outputfile = os.path.join(outputPath, wuz+'.csv')
            with open(outputfile, 'a+') as geofile:
                #geofile.write(header)
                for record in open(csvfile,'r').readlines():
                    if record.rstrip().split(',')[0] == wuz:
                        print("this record blongs to " +wuz)
                        geofile.write(record)
                    else:
                        pass

split_csvFile()


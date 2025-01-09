
# (1) First, split the occurrence file "01.geodata_cleaned_421sel.csv" into 421 seperate csv file. 
# i.e., occurrence records for each species are in a separate file；#Each file is formatted as a csv like so: 
#species,longitude,latitude
#species,longitude,latitude
#species,longitude,latitude
#species,longitude,latitude
#species,longitude,latitude
#Where there is no header line.

cd /nfs/liushuiyin/20211115_oakhyb_add/niche_reconstrc
mkdir speciesCsvFile
python 04_1.SplitCsvFile.py 01.geodata_cleaned_421sel.csv speciesCsvFile


# (2) Loop script for extracting environmental values from occurrence records
# i.e., create PNOs for each species and each environmental layer
# "extract_pointvalues_highthroughput_with_missing_data_removed_use_this.py" is from https://github.com/ryanafolk/Saxifragales_spatial_scripts/tree/master/Extract_point_values
mkdir pnos_directsampling_no_missing_data_no_point_associations
for i in ./EnvLayer/selected/*.asc; do
./extract_pointvalues_highthroughput_with_missing_data_removed_use_this.py ${i} -l ./speciesCsvFile/*.csv
done

#You should get outputs that have the environmental value in the second column and the probability of value [1/(the number of occurrences)] in the third line.

mv pnos_directsampling_no_missing_data_no_point_associations pnos

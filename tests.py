# for command line tests
import os 
import pandas as pd

#need to use 3.11.5
os.system("python3 continuous_trait_prep.py --config config.txt --output test_results_8") #initial tests

#os.system("python3 continuous_trait_prep.py --config config_v7.txt --output test_results_v7_1") #initial tests

#os.system("python3 continuous_trait_prep.py --config config_cartagene.txt --output cartagene_results_1") #make sure it works for cartagene

#os.system("python3 Script_statistics_continuous_variable.py -c cartagene/COMBINED_CATALOG_v2_9_6_JUIL2020_0.xlsx -p cartagene/gagliano_825196_data_2021_10_07_2.csv -s cartagene/samples.psam -o cartagene/gagliano_output") 

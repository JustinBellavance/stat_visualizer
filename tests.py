# for command line tests
import os 
import pandas as pd

#need to use 3.11.5

#clsa
#os.system("python3 continuous_trait_prep.py --config clsa/config_clsa.txt --output clsa/test_results") 

#os.system("python3 continuous_trait_prep.py --config config_v7.txt --output test_results_v7_1")

#CaG
#os.system("python3 continuous_trait_prep.py --config config_cartagene.txt --output cartagene_results_1") #make sure it works for cartagene

#ukbb
os.system("python3 continuous_trait_prep.py --config ukb/config_ukbb.txt --output ukb/showcase_data/test_results") 

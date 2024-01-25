import pandas as pd
import argparse
import tomllib
import numpy as np
from scipy.stats import skew, kurtosis, rankdata, norm
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import os
import re

argparser = argparse.ArgumentParser(description = 'Prepares summary stats about continuous phenotypes using a given csv file and a dictionarry.')
argparser.add_argument('-c', '--config', metavar = 'name', dest = 'config_file', type = str, required = True, help = "Config file in TOML format containing file location, sheet name and relevant column names. \n     For more information, see README")
argparser.add_argument('-o', '--output', metavar = 'name', dest = 'out_prefix', type = str, required = True, help = 'Prefix for output files.')
argparser.add_argument('-t', '--outlier_threshold', metavar = 'value', dest = 'outlier_value', type = int, default=5, help = 'Standard threshold (number of fold above and below of threshold) Fold value of the given outlier detection threshold (default 10)')

args = argparser.parse_args()


#TODO: remove --outlier threshold option
# Work to increase speed and memory usage
# Work to use data from ukbb
# Create .gitignore to make easier
# Pre-screen variables before going on to outliers?
# maybe I can run it in batches to balance between time and space, batches of 100 phenotypes? 

with open(args.config_file, mode="rb") as fp:
    TOML = tomllib.load(fp)

COLS_CODE   = TOML['COLUMN_NAMES_CODES']
COLS_MAIN   = TOML['COLUMN_NAMES_MAIN']
COLS_PHENO  = TOML['PHENOTYPE_COLS']
MISSING     = TOML['MISSING']

#these filters use the dict_categories values, but it would be better if they used the unique values from the actual data
#remove all the variables we want to ignore, if there's only 2 left, it is a binary variable
def filter_binary_variables(df, variable_missing_codes):
    label_cols = df.columns.values.tolist()
    for label in label_cols:
        df_codes = df[label].values.tolist()
        if (variable_missing_codes.get(label) != None):
            df_only_codes = np.unique(list(map(lambda x: np.nan if x in variable_missing_codes.get(label) else x, df_codes)))
        else:
            df_only_codes = np.unique(df_codes)
        df_only_codes = [x for x in df_only_codes if x == x]
        if len(df_only_codes) != 2:
            continue
        yield label
    return

#these filters use the dict_categories values, but it would be better if they used the unique values from the actual data
#remove all the variables we want to ignore, if there's only 2 left, it is a binary variable
def filter_binary_variables_2(variable_names, variable_missing_codes, contains_periods):
    index = 0
    for label in variable_names:
        index += 1
        print(f"filtering variable {index} out of {len(variable_names)}", flush = True)
        if contains_periods:
            label = f"f.{label}.0.0" 
        df = read_file(TOML["PHENOTYPE_FILE"], wanted_cols=[label])
        df_codes = df[label].values.tolist() #read in only that column
        if (variable_missing_codes.get(label) != None):
            df_only_codes = np.unique(list(map(lambda x: np.nan if x in variable_missing_codes.get(label) else x, df_codes)))
        else:
            df_only_codes = np.unique(df_codes)
        df_only_codes = [x for x in df_only_codes if x == x]
        if len(df_only_codes) != 2:
            continue
        yield label
    return

#not used anymore...
#remove all the variables we want to ignore, if there's more than 2 left, it is a categorical (or quantitative) variable
def filter_categorical_variables(df, variable_missing_codes):
    label_cols = df.columns.values.tolist()
    for label in label_cols:
        df_codes = df[label].values.tolist()
        if (variable_missing_codes.get(label) != None):
            df_only_codes = np.unique(list(map(lambda x: np.nan if x in variable_missing_codes.get(label) else x, df_codes)))
        else:
            df_only_codes = np.unique(df_codes)
        df_only_codes = [x for x in df_only_codes if x == x]
        if len(df_only_codes) > 2:
            continue
        yield label
    return

def populate_missing_codes(df):

    #get unique values in each variable name that is the same..
    dict = {'varname':[],
        'Missing codes':[]
       }

    missing = pd.DataFrame(dict)
    index = 0
    for name in df[COLS_CODE['VARIABLE']].unique():
        index += 1
        print(f"populating missing variables of column {index} out out {len(df[COLS_CODE['VARIABLE']].unique())}", flush = True)

        subset = df[df[COLS_CODE['VARIABLE']] == name]

        terms_to_avoid = '|'.join(MISSING['MISSING_CATEGORIES_SUB'])

        rows_missing_exact = subset[COLS_CODE['CATEGORY']].isin(MISSING["MISSING_CATEGORIES_EXACT"])
        
        if terms_to_avoid != "":
            rows_missing_sub = subset[COLS_CODE['CATEGORY']].str.contains(terms_to_avoid)
        else:
            rows_missing_sub = subset[COLS_CODE['CATEGORY']].str.contains(terms_to_avoid) == False

        missing_codes = subset[ rows_missing_exact | rows_missing_sub ][COLS_CODE["CODE"]].astype(str)

        codes = ",".join(missing_codes)
        missing.loc[len(missing.index)] = [name, codes] 
    return missing

def filter_continuous_variables(df_variables, variable_missing_codes, variable_names, contains_periods):

    for index, row in df_variables.iterrows():

        print(f"starting {index} out of {len(df_variables.index)}", flush = True)
        if contains_periods:
            variable = f"f.{row[COLS_MAIN['VARNAME']]}.0.0"
        else:
            variable = row[COLS_MAIN['VARNAME']]
        print(variable, flush = True)

        df_pheno_cut = read_file(TOML['PHENOTYPE_FILE'], wanted_cols =[COLS_PHENO['SEX_VAR'],variable])
        if variable not in variable_names:
            print("not in columns", flush = True)
            continue
        if (df_pheno_cut.dtypes[variable] != 'float64' and df_pheno_cut.dtypes[variable] != 'int64'):
            print("needs to be a float or int", flush = True)
            continue
        if 'RX_MED' in variable:
            print("RX_MED not relevant for GWAS", flush = True)
            continue

        unit = row[COLS_MAIN['UNITTYPE']]
        if ( COLS_MAIN['DATABASE'].strip() == "" ):
            domain = COLS_MAIN['DATABASE']
        else:
            domain = row[COLS_MAIN['DATABASE']]

        label = row[COLS_MAIN['LABEL']]

        missing_codes = variable_missing_codes.get(variable, set())

        #df_pheno_cut = df_pheno[[COLS_PHENO['SEX_VAR'], variable]].copy() # here I could just read in the column with the csv column thing :)
        
        #Removes Known missing codes
        df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: None if x in missing_codes else x )

        df_not_missing = df_pheno_cut[~df_pheno_cut.RECODED.isna()]
        if len(df_not_missing) == 0:
            continue

        mean = np.nanmean(df_pheno_cut.RECODED)
        unique_values = len(df_not_missing.RECODED.unique())
        ## Records several assessement of Outliers
        sd_details={}
        sd=np.nanstd(df_pheno_cut.RECODED)

        if sd != 0 :
            for i in range(2,int(args.outlier_value)+1):
                sd_details['[Outliers] Number of values below mean -'+str(i)+' SD'] = len(df_not_missing.RECODED[df_not_missing.RECODED < (mean-i*sd)])
                sd_details['[Outliers] Number of values above mean +'+str(i)+' SD'] = len(df_not_missing.RECODED[df_not_missing.RECODED > (mean+i*sd)])
                if sd_details['[Outliers] Number of values below mean -'+str(i)+' SD'] !=0 :
                    sd_details['[Outliers] Maximum value below mean -'+str(i)+' SD'] = max(df_not_missing.RECODED[df_not_missing.RECODED < (mean-i*sd)])
                else :
                    sd_details['[Outliers] Maximum value below mean -'+str(i)+' SD'] = 'NA'
                if sd_details['[Outliers] Number of values above mean +'+str(i)+' SD'] !=0 :
                    sd_details['[Outliers] Minimal value above mean +'+str(i)+' SD'] = min(df_not_missing.RECODED[df_not_missing.RECODED > (mean+i*sd)])
                else :
                    sd_details['[Outliers] Minimal value above mean +'+str(i)+' SD'] = 'NA'
            min_value = df_not_missing.RECODED.min()
            max_value = df_not_missing.RECODED.max()
            PERC_5,PERC95 = np.percentile(df_not_missing.RECODED , [5,95])
            S=skew(df_not_missing.RECODED, axis=0, bias=True)
        else :
            ## if SD=0 only one value is present therefor these value are irrelevant
            PERC_5='NA'
            PERC95='NA'
            S='NA'

        #If the males females are string (such as "M", "F" or "male", "female") we need to recode them to binary.
        # do this before!! not at each loop
        if (isinstance(COLS_PHENO['MALE'], str)): 
            df_not_missing.loc[:,COLS_PHENO['SEX_VAR']] = df_not_missing.loc[:,COLS_PHENO['SEX_VAR']].replace(COLS_PHENO['MALE'], 0)
        if (isinstance(COLS_PHENO['FEMALE'], str)): 
            df_not_missing.loc[:,COLS_PHENO['SEX_VAR']] = df_not_missing.loc[:,COLS_PHENO['SEX_VAR']].replace(COLS_PHENO['FEMALE'], 1)

        plot(df_not_missing, row, mean, sd) 

        # Several Statistics
        min_value = df_not_missing.RECODED.min()
        max_value = df_not_missing.RECODED.max()
        df_males = df_not_missing[df_not_missing[COLS_PHENO['SEX_VAR']] == 0]
        df_females = df_not_missing[df_not_missing[COLS_PHENO['SEX_VAR']] == 1]
        mean_value = df_not_missing.RECODED.mean()
        median_value = df_not_missing.RECODED.median()
        value_counts = df_not_missing.RECODED.value_counts()
        mode_freq = value_counts.values[0]
        mode_value = value_counts.index[0]
        n_total = len(df_not_missing)
        n_males = len(df_males)
        n_females = len(df_females)
        problem=[] ## Provides flags
        if n_total < 100 :
            problem.append('Total')
        if n_males == 0 :
            problem.append('Males')
        if n_females ==0 :
            problem.append('Females')
        if mode_freq / n_total > 0.10 : #if mode is more than 10% of values, flag problem.
            problem.append('[Statistics] Mode frequency')
        if unique_values < 30 :
            problem.append('[Statistics] N uniques values')

        print(f"done {index} out of {len(df_variables)}", flush = True)
        yield {**{
                '[Description] Variable' : variable, 
                '[Description] Domain':  domain,
                '[Description] Label' : label ,
                '[Samples] Total': n_total,
                '[Samples] Males' : n_males,
                '[Samples] Females' : n_females,
                '[Statistics] Minimum': min_value,
                '[Statistics] Maximum' : max_value, 
                '[Statistics] N uniques values': unique_values,
                '[Statistics] PERC_5' : PERC_5,
                '[Statistics] median' : median_value,
                '[Statistics] PERC_95' : PERC95,
                '[Statistics] Mean': mean_value,
                '[Statistics] Mode': mode_value,
                '[Statistics] Mode frequency': mode_freq,
                '[Hidden] Skewness': S,
                '[Statistics] SD' : sd,
                '[Hidden] problem' : problem 
            },** sd_details,**{
                '[Description] Unit': unit,
                '[Description] Missing code' : '|'.join([str(i) for i in missing_codes])
            }}

def rint(listed): # Rank Inverse Normal Transformation
    c=3/8
    rank=rankdata(listed)
    x = (rank - c) / (len(rank) - 2*c + 1)
    return norm.ppf(x)

def plot(df,row,meaned,sd): # Creates three pannel displaying the phenotype distribution
    #Raw Phenotype Distribution
    maxVal=np.nanmax(df.RECODED)
    minVal=np.nanmin(df.RECODED)
    f, axes = plt.subplots(2, 3, height_ratios=[3,1],figsize=(33, 10))

    try:
        sns.histplot(data=df, ax=axes[0,0], x="RECODED")
    except np.core._exceptions._ArrayMemoryError:
        sns.kdeplot(data=df, ax=axes[0,0], x="RECODED")

    axes[0,0].axvline(x=meaned,color='black')

    if sd !=0 :
        ulim = (maxVal + sd) if (maxVal > args.outlier_value*sd+meaned) else (args.outlier_value*sd+sd+meaned)
        llim = (minVal - sd) if (minVal < meaned-args.outlier_value*sd) else (meaned-(args.outlier_value*sd)-sd)
    else :
        ulim = maxVal + 1
        llim = minVal - 1

    axes[0,0].set_xlim(llim, ulim)
    axes[1,0].boxplot(df.RECODED,vert=False)
    axes[1,0].set_xlim(llim, ulim)
    fold=2

    ## Displays SD outlier values
    if sd != 0 :
        while (fold*sd+meaned < ulim) :
            axes[0,0].axvline(x=fold*sd+meaned,color='red')
            fold+=1
            if fold >args.outlier_value:
                break
        fold=2
        while (meaned-(fold*sd) > llim) :
            axes[0,0].axvline(x=meaned-(fold*sd),color='red')
            fold+=1
            if fold > args.outlier_value:
                break

    axes[0,0].axes.get_xaxis().set_visible(False)
    axes[1,0].axes.get_yaxis().set_visible(False)
    plt.suptitle(f'Phenotype Distribution \n {row[COLS_MAIN["LABEL"]]}',size=17)
    axes[0,0].set_ylabel('Count')
    axes[1,0].set_xlabel(row[COLS_MAIN['UNITTYPE']])

    ## Normalized Values where appropriate
    if sd !=0 :
        normalized=rint(df.loc[:,'RECODED'])
        maxValn=np.nanmax(normalized)
        minValn=np.nanmin(normalized)
        try:
            sns.histplot(data=normalized, ax=axes[0,1])
        except np.core._exceptions._ArrayMemoryError:
            sns.kdeplot(data=normalized, ax=axes[0,1])
        axes[1,1].boxplot(normalized,vert=False)
        axes[1,1].set_xlim(minValn, maxValn)
        axes[0,1].set_xlim(minValn, maxValn)
        axes[0,1].set_ylabel('Count')
        axes[1,1].set_xlabel('Normalized Unit')
        axes[0,1].axes.get_xaxis().set_visible(False)
        axes[1,1].axes.get_yaxis().set_visible(False)
        separator=(maxVal-minVal)/100 # Display purposes
        separatorN=(maxValn-minValn)/100 # Display Purposes
        ## Assess Skewness
        S=skew(df.RECODED, axis=0, bias=True)
        SS=skew(df.loc[(df.loc[:,'RECODED']>meaned-(2*sd)) & (df.loc[:,'RECODED']<meaned+(2*sd)) ,'RECODED'], axis=0, bias=True)
        SSS=skew(df.loc[(df.loc[:,'RECODED']>meaned-(3*sd)) & (df.loc[:,'RECODED']<meaned+(3*sd)) ,'RECODED'], axis=0, bias=True)
        axes[0,0].annotate("Skewness = {S:.3f} \n Skewness (without outlier 2sd) = {SS:.3f}  \n Skewness (without outlier 3sd) = {SSS:.3f} \n ".format(
        S=S,
        SS=SS,
        SSS=SSS), xy=(minVal+separator, 1), xycoords='data', xytext=(0.01, .99), textcoords='axes fraction', va='top', ha='left', fontsize='medium')
    
    ## Provides QQ plot
    sm.qqplot(df.RECODED, line ='s',ax=axes[0,2])
    plt.suptitle(f'Phenotype Distribution \n {row[COLS_MAIN["LABEL"]]}',size=13)

    #create images directory to put the plots in.
    if not os.path.exists('images'):
        os.makedirs('images')

    #save figures to that directory.
    f.savefig(f'images/{row[COLS_MAIN["VARNAME"]]}_distribution.png', dpi=300)
    plt.close(f)

# will read the file properly automatically depending on the given format (either .xlsx, .csv or .tsv) using pandas
def read_file(file_name, wanted_cols = None, only_header = False, sheet = None): 
    file_extension = os.path.splitext(file_name)[1]

    if file_extension == ".csv":
        if wanted_cols != None:
            return( pd.read_csv(file_name, header = 0, usecols=wanted_cols, low_memory = True))
        else: 
            if (only_header):
                return (pd.read_csv(file_name, index_col=0, nrows=0, low_memory = True).columns.tolist())
            else:
                return( pd.read_csv(file_name, header = 0, low_memory = True))
    elif file_extension == ".tsv" or file_extension == ".tab":
        if wanted_cols != None:
            return( pd.read_csv(file_name, header = 0, usecols=wanted_cols, low_memory = True, sep='\t'))
        else: 
            if (only_header):
                return (pd.read_csv(file_name, index_col=0, nrows=0, low_memory = True, sep='\t').columns.tolist())
            else:
                return( pd.read_csv(file_name, header = 0, low_memory = True, sep='\t'))
    elif ".xlsx":
        return(pd.DataFrame.from_dict(pd.read_excel(file_name, sheet_name = sheet, header = 0)))
    
    return None

#because ukb phenotype file is formatted as f.<fieldID>.<visit#>.<repeatedmeasurement#> 
#We need to change them to only get the <fieldID> or else it won't correspond the with dictionnary.
#this script won't affect any other phenotype files because they don't have periods in the header.
def clean_ukb_format(headers):

    contains_periods = "." in headers[0]
    #if the variable names contain periods (like with ukb format) take the variable of interest
    if (contains_periods):
        headers_list = []
        for header in headers:
            if header.split(".")[2] == "0" and header.split(".")[3] == "0":
                headers_list.append(header.split(".")[1])
        headers_cleaned = set(headers_list)
    else:
        headers_cleaned = set(headers)

    # remove repetitive var names
    return headers_cleaned, contains_periods

if __name__ == '__main__':

    #if theres's seperate coding file (such as with ukb), read on only that full file
    print("Reading code files", flush = True)
    if TOML['CODING_FILE'] != "":
        codes = read_file(TOML['CODING_FILE'])
    else:
        codes = read_file(TOML['DICT_FILE'], sheet = TOML['SHEET_NAME_CODES'])

    #remove columsn that might be empty or missing a header
    print("removing missing values in code columns")
    subset_values = []
    for value in COLS_CODE:
        if (COLS_CODE[value] != ""):
            subset_values.append(COLS_CODE[value])

    print(TOML['SHEET_NAME_CODES'])
    codes = codes[subset_values]

    #read in phenotypes    
    print("Reading in phenotype header", flush = True)
    variable_names = read_file(TOML['PHENOTYPE_FILE'], only_header = True)
    print(variable_names, flush = True)
    variable_names, contains_periods = clean_ukb_format(variable_names)
    print(variable_names, flush = True)

    #make sure all the code labels that we will be using are strings.
    print("Converting code labels to string", flush = True)
    codes[COLS_CODE['CATEGORY']] = codes[COLS_CODE['CATEGORY']].astype("string")
    print(codes[COLS_CODE['CATEGORY']].head(), flush = True)

    #run function to create dataframe in needed format. E.g.:
    #   varname Missing codes
    # 0 ACUTE_RENAL_FAIL_ONSET_AGE  -9, -7, 77, 88, 99
    # 1 ADDICTION_DISORDER_ONSET_AGE    -9, -7, 77, 88, 99 
    print("Populating missing codes... ", flush = True)
    df_missing_codes = populate_missing_codes(codes)
    print("Populating missing codes... Done.", flush = True)

    #remove empty missing codes (as they might appear in code files)
    print(df_missing_codes)
    df_missing_codes = df_missing_codes[df_missing_codes.iloc[:,1] != ""]
    df_missing_codes.to_csv(f'{args.out_prefix}.missing_codes.csv')
    
    #remove all non digits, make all digits floats.
    #this will make a json that has all the missing codes for a given variable
    variable_missing_codes = {row['varname']: set(None if re.search('[^0-9]+', x) else float(x) for x in row['Missing codes'].split(',')) for index, row in df_missing_codes.iterrows()}

    #binary_variables = set(filter_binary_variables(df_pheno, variable_missing_codes)) 
    print("binary variable filtering ...", flush = True)
    binary_variables_2 = set(filter_binary_variables_2(variable_names, variable_missing_codes, contains_periods))
    print("binary variable filtering ... done", flush = True)
    #categorical_variables = set(filter_categorical_variables(df_pheno, variable_missing_codes))

    trait_dict = read_file(TOML['DICT_FILE'], sheet = TOML['SHEET_NAME_MAIN'])

    trait_dict[COLS_MAIN['UNITTYPE']] = trait_dict[COLS_MAIN['UNITTYPE']].fillna('') #if there is a missing space, treat as NA

    #remove binary (to keep only quantitative variables)
    trait_dict = trait_dict[~trait_dict[COLS_MAIN['VARNAME']].isin(binary_variables_2)]
    #trait_dict = trait_dict[~trait_dict[COLS_MAIN['VARNAME']].isin(categorical_variables)]
        
    #populate variables for final dataframe.
    print("Starting filter_continous variables...", flush = True)
    #is there a way to write this line by line as a file, and do the .to_json in more memory efficient way on the command line or something? (probably)
    continuous_variables_final = pd.DataFrame(filter_continuous_variables(trait_dict, variable_missing_codes, variable_names, contains_periods)) 
    print("Done filtering continuous variables, saving to JSON..", flush = True)

    # filter_continuous_variables(trait_dict, variable_missing_codes, variable_names, contains_periods).to_json(f'{args.out_prefix}.summary.json', orient='records', lines = True, mode = 'a')

    continuous_variables_final.to_json(f'{args.out_prefix}.summary.json', orient='records', lines = True, mode = 'a')
    print("script done!")

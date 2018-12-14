import pandas as pd
import numpy as np
import regex as re
import warnings
import glob
import os
from functools import reduce

def geno(glob_path, clean, chars):
    """ Processes ARIBA data
    """

    #initialise empty list for assembled,ref_seq dataframes
    df_list = []
    dflist1 = []

    #initialise empty list to append dfs
    frame = []    

    #initialise empty list for unique names from all tables
    names = []

    #read in csv files as per defined user path and print a message listing each file found
    for n, csv in enumerate(glob.glob('{}/*.csv'.format(glob_path))):

        print('\t Found csv file {}: {}...'.format(n+1, csv))

        df = pd.read_csv(csv)

        #check columns to see if the file is in the correct format, giving an error message if not
        assembled = (re.search(r".*\.ref_seq", df.columns[1]))
        ref = (re.search(r".*\.ref_seq", df.columns[2]))

        correctcols = [assembled,ref]

        if any(correctcols) == False:
            print("\nColumn names in the above file appear to be incorrect\n\nDid you generate the ARIBA summaries as follows?\n\tariba summary --cluster_cols assembled,ref_seq <summary_prefix> <ARIBA_output/report.tsv>")
            print("\nAlternatively you may have saved your MLST file as a CSV rather than a .txt or have other CSV files in this directory other than files to be processed by ARIBAlord...\n")
            exit()

        #clean sample name
        if clean:
            df['name'] = df['name'].replace("{}.*".format(chars), "", regex = True)
        else:
            df['name'] = df['name'].replace("_.*", "", regex = True)

        #change index to name column
        df.set_index('name',inplace=True)

        #strip whitespace from column heads
        df_clean = df.rename(columns=lambda x: x.strip())

        #check read csv files to see if they are phylogroup data - this is processed differently
        phylogroup = [(re.search(r"tspE4|arpA|chuA|yjaA", i)) for i in df_clean]

        #convert string values to representative integers. Values are different for phylogroup
        if any(phylogroup):
            cats = ['yes', 'no', 'yes_nonunique', 'partial', 'interrupted', 'fragmented']
            values = [1, 0, 1, 1, 1, 1]
            df_test = df_clean.replace(to_replace = cats, value = values)
        else:
            cats = ['yes', 'no', 'yes_nonunique', 'partial', 'interrupted', 'fragmented']
            values = [1, 0, 1, 0, 0, 0]
            df_test = df_clean.replace(to_replace = cats, value = values)

        #initialise empty list for assembled,ref_seq dataframes
        df_list = []

        #define number of columns, convert to a list, iterate this list and bundle into batches of two
        new_df_length = int(len(df.columns))
        split = list(range(new_df_length))
        it = iter(split)
        zip_it = zip(it, it)

        #for the two columns associated with each gene hit
        #save these columns and their index to a data frame
        #then append this dataframe to a list of dataframes
        for i in zip_it:
            df_pair = df_test.iloc[:,[i[0],i[1]]]
            df_list.append(df_pair)
           

        #initialise empty list
        frame = []   

        #for dataframe in df_list
        #pivot the tables based so gene hits are colnames and 'assembled' status are observations
        #remove trash column that is generated in column 0
        #change NaNs to zeroes
        #append these modified dataframes to new list of dataframes called frame
        for i in df_list :
            data_pivot = (i.pivot(values = i.columns[0], columns = i.columns[1]))
            if len(data_pivot.columns) == 1:
                data_pivot['NaN'] = np.nan
                data_pivot = data_pivot.iloc[:,[1,0]]
            data_pivot = data_pivot.iloc[:,1:]
            data_pivot = data_pivot.fillna(0)
            frame.append(data_pivot)
            

                       

        #make list of original sample names using index of original df    
        base = df.iloc[:,0:0]

        #join new dataframes from frames to base dataframe
        data_all = base.join(frame)

        dflist1.append(data_all)

        names.append(base)

    return (dflist1, names)

def simple_clean(dflist):
#regular expressions to clean column names
    dflist2 = []
    
    for df in dflist:

        #patterns used to recognise spreadsheets as being of serotype, custom, virulence etc
        #to apply the correct regexs to clean column names
        phylogroup = [(re.search(r"^arpA|^tspE4\.C2|^yjaA", i)) for i in df]
        virulence = [(re.search(r"^gad", i)) for i in df]
        resistance = [(re.search(r"^dfrA|^tet|^bla|^sul|^aadA", i)) for i in df]
        plasmid = [(re.search(r"^Inc", i)) for i in df]
        insertion = [(re.search(r"^IS", i)) for i in df]
        custom = [(re.search(r"^fimH", i)) for i in df]
        serotype = [(re.search(r"^wzy|^wzx|^wzt|^wzm|^fliC", i)) for i in df]
        srst2serotype = [(re.search(r"^[0-9].*wzy|^[0-9].*wzx|^[0-9].*wzt|^[0-9].*wzm|^[0-9].*fliC", i)) for i in df]
        MLST = [(re.search(r"^adk|^fumC|^gyrB", i)) for i in df]
        
        if any(phylogroup):
            df = df.rename(columns=lambda x: re.sub('(^arpA|^tspE4\.C2|^yjaA|^chuA)',r'phylogroup_\1',x))

        elif any(custom):
            df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
            df = df.rename(columns=lambda x: re.sub('^(.*)(NC_)?_.*',r'v_\1',x))
            df = df.rename(columns=lambda x: re.sub('ipaH','v_ipaH',x))
            df = df.rename(columns=lambda x: re.sub('^v_IS','i_IS',x))
            df = df.rename(columns=lambda x: re.sub('IS1_A','IS1',x))
            df = df.rename(columns=lambda x: re.sub('^v_int','i_int',x))
            df = df.rename(columns=lambda x: re.sub('^v_(mer|ter|sil|pco|czc|scs)',r'r_\1',x))
            df = df.rename(columns=lambda x: re.sub('_$','',x))
            df = df.rename(columns=lambda x: re.sub('_NC$','',x))
            df = df.rename(columns=lambda x: re.sub('v_kpsMT_II__K2_CP000468.1__APEC','v_kpsMT_II',x))
            df = df.rename(columns=lambda x: re.sub('v_kpsMT_II__K1_AE014075.1','v_kpsMT_II',x))
            df = df.rename(columns=lambda x: re.sub('v_kpsMT_III__K54','v_kpsMT_III',x))
            df = df.rename(columns=lambda x: re.sub('v_kpsMT_II__K5_CP002212.1_clone_D','v_kpsMT_II',x))
            df = df.rename(columns=lambda x: re.sub('afaE','v_afaE',x))

        elif any(virulence):
            df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
            df = df.rename(columns=lambda x: re.sub('^(.*)',r'v_\1',x))

        elif any(serotype):
            print('\n serotypefinder-serotype detected\n')
            df = df.rename(columns=lambda x: re.sub('(^[^_]+)_[^_]+_[^_]+_(.*)',r'sero_\1:\2',x))
            

        elif any(srst2serotype):
            df = df.rename(columns=lambda x: re.sub('[0-9]__[A-Za-z]+__([^-]+)_([^_]+)__.*',r'sero_\1:\2',x))
            print('\n srst2 serotype detected\n')

        elif any(resistance):
            df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
            df = df.rename(columns=lambda x: re.sub('^(.*)',r'r_\1',x))
            df = df.rename(columns=lambda x: re.sub('aph_3____Ib','strA',x))
            df = df.rename(columns=lambda x: re.sub('aph_6__Id','strB',x))
            df = df.rename(columns=lambda x: re.sub('aac_3__IId','aac(3)-IId',x))
            df = df.rename(columns=lambda x: re.sub('r_aph_3___Ia','r_aph(3\')Ia',x))
            df = df.rename(columns=lambda x: re.sub('r_aph_4__Ia','r_aph(4\')Ia',x))
            df = df.rename(columns=lambda x: re.sub('blaTEM_','blaTEM-',x))
            df = df.rename(columns=lambda x: re.sub('lnu_([A-Z])_',r'lnu(\1)',x))
            df = df.rename(columns=lambda x: re.sub('mph_([A-Z])_',r'mph(\1)',x))
            df = df.rename(columns=lambda x: re.sub('tet_([A-Z])_',r'tet(\1)',x))


        elif any(plasmid):
            df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
            df = df.rename(columns=lambda x: re.sub('^FIA',r'IncFIA',x))
            df = df.rename(columns=lambda x: re.sub('^(.*)',r'p_\1',x))
            df = df.rename(columns=lambda x: re.sub('(p_IncF[A-Z]+)_.*',r'\1',x))
            df = df.rename(columns=lambda x: re.sub('(p_Col_MG828)_',r'\1',x))
            df = df.rename(columns=lambda x: re.sub('p_IncB_O_K_Z','p_IncB/O/K/Z',x))
            df = df.rename(columns=lambda x: re.sub('p_IncHI1B_CIT_','p_IncHI1B',x))

        elif any(insertion):
            df = df.rename(columns=lambda x: re.sub('([^-])_.*',r'\1',x))
            df = df.rename(columns=lambda x: re.sub('^(IS.*)',r'i_\1',x))

        elif any(MLST):
            warnings.warn('CSV file appears to contain MLST alleles.', Warning)
        else:
            warnings.warn('CSV file contents could not be identified and therefore may have sloppy column headers.\n This may be due to the use of a custom database', Warning)
        
        df = df.rename(columns=lambda x: re.sub('.*name','name',x))

        df = df.groupby(df.columns, axis=1).sum()

        dflist2.append(df)
    
    return dflist2


def mlst(glob_path):

    for n, txt in enumerate(glob.glob('{}/*.tsv'.format(glob_path))):

        print('\tFound MLST file: {}...'.format(txt))

        MLST = pd.read_csv(txt, delimiter='\t')

        MLST_pattern = [(re.search(r"^adk|^fumC|^gyrB|^icd|^mdh|purA|recA", i)) for i in MLST]

        if any(MLST_pattern) :
            MLST = MLST.rename(columns=lambda x: re.sub(r'(adk|^fumC|^gyrB|^icd|^mdh|^purA|^recA)',r'MLST_\1',x))
        else :
            warnings.warn('MLST text file does not contain expected values...', Warning)

        MLST.columns.values[-1] = 'name'

        print(MLST)

        if args.clean:
            MLST['name'] = MLST['name'].replace("{}.*".format(chars), "", regex = True)
        else:
            MLST['name'] = MLST['name'].replace("_.*", "", regex = True)

        MLST = MLST[~MLST["ST"].str.contains("ST")]

        print(MLST.head())

        mlst_table = MLST

        return mlst_table
    

def sero(table, simple_csv=True):
    
    pd.options.mode.chained_assignment = None

    table = table.reset_index()

    table = table.filter(regex=r'(^name|^sero_)', axis=1)


    table = table.rename(columns=lambda x: re.sub(r'^sero_','',x))
    
   
    #melt df
    table = pd.melt(table, id_vars='name')

    #filter non 1s
    table = table[table.value == 1]

    #split string of fliC_H4, for example, to two columns, one with fliC and one with H4
    table['simple'] = table.variable.str.split(':').str[0]
    table['type'] = table.variable.str.split(':').str[1]
  
    #remove unwanted columns
    table = table.iloc[:,[0,3,4]]

    #group by name and gene type ('simple') and concatenate together multiple hits for the same gene type ('simple')
    table['O_or_H'] = table.groupby(['name','simple'])['type'].transform(lambda x: '-'.join(x))

    #remove unwanted columns
    table = table.iloc[:,[0,1,3]]

    # #drop duplicate rows
    table.drop_duplicates(inplace=True)

    # #pivot table
    table = table.pivot(index = 'name', columns = 'simple', values = 'O_or_H')


    table['wzy'] = table['wzy'].replace("^[0-9]+_", "", regex = True)

    no_wzm_wzt = [(re.search(r"wzt|wzm", i)) for i in table]
    no_wzx_wzy = [(re.search(r"wzx|wzy", i)) for i in table]

    if not any(no_wzm_wzt):
        table['wzm'] = np.nan
        table['wzt'] = np.nan

    if not any(no_wzx_wzy):
        table['wzx'] = np.nan
        table['wzy'] = np.nan

    #create empty columns with NaNs - i believe this is redundant...
    table['O1cat'] = np.nan
    table['O2cat'] = np.nan
    table['O_type'] = np.nan
    table['H_type'] = np.nan

    non_fliCs_present = [(re.search(r"^(flnA|flmA|fllA|flkA)", i)) for i in table]

    #if no non-fliC H hits then set H_type to fliC, otherwise set to fliC and add an asterisk
    if non_fliCs_present == False:
        table['H_type'] = table['fliC']+"*"
    else:
        table['H_type'] = table['fliC']

    #replace null with blank
    table = table.where((pd.notnull(table)), str('*'))

    #add two new columns combining the gene hits for wzm/wzt and wzx/wzy
    table['O1cat'] = table['wzm'].map(str) + "*" + '/' + table['wzt'] + "*"
    table['O2cat'] = table['wzx'].map(str) + "*" + '/' + table['wzy'] + "*"

    #There are some tricky combinations of wzx/wzy that choke up our serotyping processing.
    #Below are regular expressions that process such cases
    table['O2cat'] = table['O2cat'].replace("O17_O77\*\/O17_O44\*","O17/O44/O77", regex = True)
    table['O2cat'] = table['O2cat'].replace("O50_O2\*\/O2\*","O2", regex = True)
    table['O2cat'] = table['O2cat'].replace("O50_O2\*\/O50\*","O50", regex = True)
    table['O2cat'] = table['O2cat'].replace("O118_O151\*\/O118\*","O118", regex = True)
    table['O2cat'] = table['O2cat'].replace("O118_O151\*\/O151\*","O151", regex = True)
    table['O2cat'] = table['O2cat'].replace("O135\*\/O13_O135\*","O135", regex = True)
    table['O2cat'] = table['O2cat'].replace("O13\*\/O13_O135\*","O13", regex = True)
    table['O2cat'] = table['O2cat'].replace("O123\*\/O123_O186\*","O123", regex = True)
    table['O2cat'] = table['O2cat'].replace("O186\*\/O123_O186\*","O186", regex = True)
    table['O2cat'] = table['O2cat'].replace("O169\*\/O169_O183\*","O169", regex = True)
    table['O2cat'] = table['O2cat'].replace("O183\*\/O169_O183\*","O183", regex = True)
    table['O2cat'] = table['O2cat'].replace("O141ab_O141ac\*\/O141ab\*","O141b", regex = True)
    table['O2cat'] = table['O2cat'].replace("O141ab_O141ac\*\/O141ac\*","O141c", regex = True)

    #if wzm==wzt then set O1 (O_type 1) to wzm, otherwise set O1 to O1cat
    table['O1'] = np.where(table['wzm']==table['wzt'], table['wzm'], table['O1cat'])

    #if wzx==wzy then set O2 (O_type 2) to wzx, otherwise set O2 to O2cat
    table['O2'] = np.where(table['wzx']==table['wzy'], table['wzx'], table['O2cat'])

    #replace */* with an empty cell
    table = table.replace("\*/\*$", "", regex = True)

    #combined O1 and O2 with a slash between them
    table['O_type'] = table['O1'].map(str) + '/' + table['O2']

    #replace separating strings that exist in absence of gene hits
    table = table.replace("\*\*\/\*\*", "", regex = True)
    table = table.replace("\*\/\*", "", regex = True)
    table = table.replace("^\*\/", "", regex = True)
    table = table.replace("\/\*", "", regex = True)


    #replace non-hits with ONT
    table['O_type'] = table['O_type'].replace("^$", "ONT", regex = True)
    table['H_type'] = table['H_type'].replace("^\*?$", "HNT", regex = True)

    #combine O and H type into O:H format in column OH_type
    table['OH_type'] = table['O_type'].map(str) + ':' + table['H_type']

    table = table.reset_index()

    #make simple table of EcOH data
    EcOH = table.loc[:,['name','O_type','H_type', 'OH_type']]

    #create a more informative table
    seroframe = table
    
    if simple_csv:
        return EcOH
    else:
        return table


def phylog(table):
    
    #Pull out columns with names 'name' or starting with 'phylogroup'
    phylogroup = table.filter(regex=r'(^name|^phylogroup_)', axis=1)

    arpA_pres = [(re.search(r"^phylogroup_arpA", i)) for i in phylogroup]
    tspE4_pres = [(re.search(r"^phylogroup_tspE4\.C2", i)) for i in phylogroup]
    yjaA_pres = [(re.search(r"^phylogroup_yjaA", i)) for i in phylogroup]
    chuA_pres = [(re.search(r"^phylogroup_chuA", i)) for i in phylogroup]

    if any(arpA_pres):
        print("arpA detected")
    else:
        print("arpA not detected, column of zeros added for processing...")
        phylogroup.loc[:,'phylogroup_arpA'] = 0.0

    if any(tspE4_pres):
        print("tspE4.C2 detected")
    else:
        print("tspE4.C2 not detected, column of zeros added for processing...")
        phylogroup.loc[:,'phylogroup_tspE4.C2'] = 0.0

    if any(yjaA_pres):
        print("yjaA detected")
    else:
        print("yjaA not detected, column of zeros added for processing...")
        phylogroup.loc[:,'phylogroup_yjaA'] = 0.0

    if any(chuA_pres):
        print("chuA detected")
    else:
        print("chuA not detected, column of zeros added for processing...")
        phylogroup.loc[:,'phylogroup_chuA'] = 0.0
    

    #Trim 'phylogroup' prefix
    phylogroup = phylogroup.rename(columns=lambda x: re.sub(r'^phylogroup_','',x))

    #Set index to name to simplify processing
    phylogroup.set_index('name')

    #Subset samples based on gene carriage as per phylogroup classification
    B2_or_D = phylogroup.loc[phylogroup['chuA'] == 1]
    A_or_B1 = phylogroup.loc[phylogroup['chuA'] == 0]
    B2 = B2_or_D.loc[phylogroup['yjaA'] == 1]
    D = B2_or_D.loc[phylogroup['yjaA'] == 0]
    B1 = A_or_B1.loc[phylogroup['tspE4.C2'] == 1]
    A = A_or_B1.loc[phylogroup['tspE4.C2'] == 0]



    #Initialise empty list and append the subsetted dataframes to it
    listdfs = []
    listdfs = listdfs.append([A,B1,B2,D])

    #Classify samples based on the dataframe they have been subset into
    A['phylogroup'] = 'A'
    B1['phylogroup'] = 'B1'
    B2['phylogroup'] = 'B2'
    D['phylogroup'] = 'D'

    #Combine dataframe and gene columns
    combined = pd.DataFrame()
    combined = pd.DataFrame(combined.append([A,B1,B2,D]))
    summary = combined.loc[:,['name','phylogroup']]

    return summary

def report(simple):
    #Generate counts of STs and Novel's detected by ARIBA
    STcount = simple['ST'].replace('\*','',regex=True)
    Novelcount = STcount[STcount == 'Novel']
    STcount = STcount[STcount != 'Novel']
    STcount = len(STcount.unique())
    Novelcount = len(Novelcount)

    #Calculate genenum from columns in dataframe and samplenum from number of rows in dataframe
    genenum = simple.shape[1]
    samplenum = simple.shape[0]

    #Sort columns by their names
    simple = simple.reindex(sorted(simple.columns), axis=1)

    #Remove columns with sum zero
    simple = simple.loc[:, (simple != 0).any(axis=0)]

    #Filter virulence columns
    simplevir = simple.filter(regex="^v_")
    #Convert anything greater than a one to a one to facilitate sum
    simplevir[simplevir > 1] = 1
    #filter rows based on which the highest count of virulence genes
    simplevirsum = simplevir[simplevir.sum(axis=1) == simplevir.sum(axis=1).max()]

    #As above for resistance
    simpleres = simple.filter(regex="^r_")
    simpleres[simpleres > 1] = 1
    simpleressum = simpleres[simpleres.sum(axis=1) == simpleres.sum(axis=1).max()]

    #As above for plasmid
    simpleplas = simple.filter(regex=r"^(p_Inc|p_p0111)")
    simpleplas[simpleplas > 1] = 1
    simpleplassum = simpleplas[simpleplas.sum(axis=1) == simpleplas.sum(axis=1).max()]

    #As above for insertion/integrase
    simpleins = simple.filter(regex="^i_")
    simpleins[simpleins > 1] = 1
    simpleinssum = simpleins[simpleins.sum(axis=1) == simpleins.sum(axis=1).max()]

    #Initialises an empty list we can append to facilitate report writing
    towrite = []

    #Generating statistics to report on
    towrite.append('{} read-sets analysed'.format(samplenum))
    towrite.append('\n{} sequence types identified'.format(STcount))
    towrite.append('\n{} samples may constitute novel STs'.format(Novelcount))
    towrite.append('\n{} genes identified'.format(genenum))

    towrite.append('\nResistance genes detected include:')
    for i in simpleres.columns:
        toappend = str(re.sub('r_','',i) + "\t" + str(int(sum(simple.loc[:,i]))))
        towrite.append(toappend)

    #A series of for loops that print the sample phylonames for samples that have max counts of genes of various associations
    #Vir Max
    towrite.append("\n Max virulence count ("+str(simplevir.sum(axis=1).max()) + ") detected in the following samples: ")
    simplevirsum.reset_index(inplace=True)
    for i in simplevirsum['phyloname']:
        towrite.append(i)
    #Res Max
    towrite.append("\n Max resistance count ("+str(simpleres.sum(axis=1).max()) + ") detected in the following samples: ")
    simpleressum.reset_index(inplace=True)
    for i in simpleressum['phyloname']:
        towrite.append(i)
    #Plas Max
    towrite.append("\n Max plasmid count ("+str(simpleplas.sum(axis=1).max()) + ") detected in the following samples: ")
    simpleplassum.reset_index(inplace=True)
    for i in simpleplassum['phyloname']:
        towrite.append(i)
    #Insertion/Integrase Max
    towrite.append("\n Max insertion count ("+str(simpleins.sum(axis=1).max()) + ") detected in the following samples: ")
    simpleinssum.reset_index(inplace=True)
    for i in simpleinssum['phyloname']:
        towrite.append(i)

    #Open a file for us to write our report
    report = open('{}_report.txt'.format(args.output_file), 'w')

    #Write items from towrite to report file
    for item in towrite:
        report.write("{}\n".format(item))

    #Close file and print outfile name
    report.close()
    print("\nARIBAlord report written to {}_report.txt".format(args.output_file))


if __name__ == '__main__':
    import argparse

parser = argparse.ArgumentParser(description='Process and join all genotypic and phylogenetic ARIBA output')
parser.add_argument('input_dir', help='Input directory containing csv files')
parser.add_argument('output_file', help='Output filename')
parser.add_argument('--clean', action='store_true', help='Use to clean anything after _R1 from filenames')
parser.add_argument('-chars','--characters', action='store_true', help='Pattern to trim after')
args = parser.parse_args()

if args.clean:
    print('\tRunning geno on CSV files in {} \nTrimming characters in column \'name\' after \'{}\'...'.format(args.input_dir, args.characters))
    dflist1, names = geno(args.input_dir, args.clean, args.characters)
else:
    print('Running geno on CSV files in {}... \nArgument \'--trim\' not used.'.format(args.input_dir))
    dflist1, names = geno(args.input_dir, args.clean, args.characters)

#Used to generate unique list of names from across all CSVs
catnames = pd.concat(names)
catnames.reset_index(inplace=True)
names = catnames.drop_duplicates(keep='first')

print('\nCleaning headers of processed CSV files, adding identifiers prefixes')
dflist2 = simple_clean(dflist1)

#Check for presence of MLST, EcOH and Phylogroup files

if glob.glob('{}/*MLST*.tsv'.format(args.input_dir)):
    MLST_present = True
else:
    MLST_present = False

if glob.glob('{}/*EcOH*.csv'.format(args.input_dir)):
    EcOH_present = True
else:
    EcOH_present = False

if glob.glob('{}/*hylogroup*.csv'.format(args.input_dir)):
    Phylogroup_present = True
else:
    Phylogroup_present = False

Phylogenetic_presents = [MLST_present, EcOH_present, Phylogroup_present]

if MLST_present:
    print('\nGenerating MLST table:')
    mlst_table = mlst(args.input_dir)

    #Set index to name and take only name and ST from MLST table
    mlst_table = mlst_table.set_index('name')
    mlst_simple = mlst_table.iloc[:,0].to_frame()

    #Join the processed tables generated in geno
    full = mlst_table.join(dflist1)
    simple = mlst_simple.join(dflist2)

else:
    print('\nNo MLST table detected; skipping MLST...')
    names = names.set_index('name')
    full = names.join(dflist1)
    simple = names.join(dflist2)
    full = full.reset_index()
    simple = simple.reset_index()

if EcOH_present:
    print('\nEcOH detected')
    EcOH = sero(simple, simple_csv=True)
    serotable = sero(simple, simple_csv=False)
else:
    print('\nNo EcOH table detected; skipping EcOH...')

if Phylogroup_present:
    print('\nPhylogroup detected')
    #Process phylogroup data
    simple = simple.reset_index()
    phylogroup = phylog(simple)
else:
    print('\nNo Phylogroup table detected; skipping Phylogroup...')

    #Merge gene hits with serotype and phylogroup data
if EcOH_present:
    simple = pd.merge(simple,EcOH,how='outer',on='name')
if Phylogroup_present:
    simple = pd.merge(simple,phylogroup,how='outer',on='name')

#Remove columns that contain raw sero/phylogroup hits
if EcOH_present:
    simple = simple.loc[:, ~simple.columns.str.startswith('sero_')]
if Phylogroup_present:
    simple = simple.loc[:, ~simple.columns.str.startswith('phylogroup_')]


if EcOH_present:
    simple['O_type'] = simple['O_type'].fillna(value = 'ONT')
    simple['H_type'] = simple['H_type'].fillna(value = 'HNT')
    simple['OH_type'] = simple['OH_type'].fillna(value = 'ONT:HNT')

simple = simple.fillna(value = 0)

#Create column 'phyloname' from name and phylogenetic data
# if Phylogenetic_presents == [False, False, False]:
#     simple['phyloname'] = simple['name']
#     simple = simple.set_index(['name','phyloname'])
# if Phylogenetic_presents == [False, False, True]:
#     simple['phyloname'] = simple['name']+simple['phylogroup']
#     simple = simple.set_index(['name','phyloname','phylogroup'])
# if Phylogenetic_presents == [False, True, False]:
#     simple['phyloname'] = simple['name']+"-"+simple['OH_type']
#     simple = simple.set_index(['name','phyloname','O_type','H_type','OH_type'])
# if Phylogenetic_presents == [False, True, True]:
#     simple['phyloname'] = simple['name']+"-"+simple['phylogroup']+"-"+simple['OH_type']
#     simple = simple.set_index(['name','phyloname','phylogroup','O_type','H_type','OH_type'])
# if Phylogenetic_presents == [True, False, True]:
#     simple['phyloname'] = simple['name']+"-ST"+simple['ST']+"-"+simple['OH_type']
#     simple = simple.set_index(['name','phyloname','ST','O_type','H_type','OH_type'])
# if Phylogenetic_presents == [True, True, True]:
#     simple['phyloname'] = simple['name']+"-ST"+simple['ST']+"-"+simple['phylogroup']+"-"+simple['OH_type']
#     simple = simple.set_index(['name','phyloname','ST','phylogroup','O_type','H_type','OH_type'])

simple = simple.reset_index()

index_present = [(re.search(r"^index$", i)) for i in simple]

if any(index_present):
    simple = simple.drop(labels = 'index', axis = 1)

simple = simple.loc[:, (simple != 0).any(axis=0)]

#Write final tables to CSV
print('\nWriting simplified ARIBA table to {}_simple.csv'.format(args.output_file))
simple.to_csv(args.output_file+'_simple.csv', index = False)
print('\nWriting full ARIBA table to {}_full.csv\n'.format(args.output_file))
full.to_csv(args.output_file+'_full.csv', index = False)
if EcOH_present:
    print('\nWriting full serotype table {}_sero.csv\n'.format(args.output_file))
    serotable.to_csv(args.output_file+'_sero.csv', index = False)

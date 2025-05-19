import pandas as pd
import numpy as np
import re
from dotenv import load_dotenv
import os


#RH

def process_RH(df,subset='BRCA1_BRCA2_set'):

    # Combine the first three rows into one header
    new_header = df.iloc[:3].apply(lambda x: '_'.join(x), axis=0)
    df.columns = new_header
    df = df.drop(index=[0, 1, 2]) # Drop the first three rows since they're now part of the header

    df.reset_index(drop=True, inplace=True)
    df.rename(columns={'Cell line_Condition_sgRNA sequence': 'Guide'}, inplace=True)

    bes= pd.read_excel("data/Base Editing Screens Samplesheet.xlsx",sheet_name='BE sgRNA library_RH' )
    if subset in bes['Subset'].unique():
        bes= bes[bes['Subset'] == subset]       
    else:
        raise print('Subset Doesnt Exist')
        
    
    df2_selected = bes[['sgRNA sequence', 'Editor','GeneID']]
    df2_selected.rename(columns={'sgRNA sequence': 'Guide'}, inplace=True)
    df2_selected.rename(columns={'GeneID': 'Gene'}, inplace=True)
    df2_selected['Editor'] = df2_selected['Editor'].replace({'A-G': 'ABE', 'C-T': 'CBE'})

    assert bes['sgRNA sequence'].nunique()==bes.shape[0]
    table=df.merge(df2_selected, on="Guide",how="left").drop_duplicates() #here we can merge on guide only because for each subset there are only unique guides
    new_order = ['Guide', 'Gene', 'Editor'] + [col for col in table.columns if col not in ["Guide", "Gene", "Editor"]]
    df = table[new_order]
    df= df.sort_values(by="Gene")
    df = df.reset_index(drop=True)
    df.rename(columns={'pDNA_pDNA_pDNA': 'pDNA'}, inplace=True)
    return df


#GR

def process_GR(df):
    for col in df.columns:
        if 'pDNA' in col:
            df = df.rename(columns={col: 'pDNA'}, inplace=False)

    df = df.rename(columns={'Construct Barcode': 'Guide'}, inplace=False)

    
    editor_types = set()
    data_cols = []
    for col in df.columns:
        match_old = re.match(r'([A-Za-z0-9]+)_([A-Za-z]+)_([A-Za-z]+)', col)
        match_new = re.search(r'_([A-Za-z]+)_', col)

        if match_old:
            editor_types.add(match_old.group(2))
            data_cols.append(col)
        elif match_new:
            editor_types.add(match_new.group(1))
            data_cols.append(col)

    if not editor_types:
        return df #Return the original dataframe if no editor types were found.

    reshaped_dfs = []
    for editor in editor_types:
        editor_cols = [col for col in data_cols if f"_{editor}_" in col or f"{editor}_" in col]
        if editor_cols:
            editor_df = df[['Guide', 'pDNA'] + editor_cols].copy()
            editor_df['Editor'] = editor
            if f"_{editor}_" in editor_cols[0]:
                editor_df.columns = ['Guide', 'pDNA'] + [re.sub(f"_{editor}_", "", col) for col in editor_cols] + ['Editor']
            else:
                editor_df.columns = ['Guide', 'pDNA'] + [re.sub(f"{editor}_", "", col) for col in editor_cols] + ['Editor']
            reshaped_dfs.append(editor_df)

    merged_df = pd.concat(reshaped_dfs, ignore_index=True)

    rep_columns = []
    for column in merged_df.columns:
        if re.search(r'Rep[A-Za-z]', column) or re.search(r'rep[A-Za-z]', column):
            rep_columns.append(column)

    rep_columns_with_underscore = [re.sub(r'([Rr]ep[A-Za-z])', r'_\1', col) for col in rep_columns]

    reordered_columns = ['Guide','Editor'] + rep_columns_with_underscore + ['pDNA']

    # Check if all reordered columns are in the merged_df
    if all(col in merged_df.columns for col in reordered_columns):
        merged_df = merged_df[reordered_columns]
    else:
        print("Warning: Some reordered columns not found in DataFrame.")

    merged_df = merged_df.replace(0, np.nan)
    merged_df = merged_df.sort_values(by=['Guide','Editor'])

    new_columns = []
    for col in merged_df.columns:
        if re.search(r'[Rr]ep[A-Za-z]', col):
            new_col = re.sub(r'([Rr]ep[A-Za-z])', r'_\1', col)
        else:
            new_col = col
        new_columns.append(new_col)
    merged_df.columns = new_columns


    ## HERE WE GET THE GENES FROM BES
    bes= pd.read_excel("data/Base Editing Screens Samplesheet.xlsx",sheet_name='BE sgRNA library_GR' )
    df2_selected = bes[['sgRNA sequence','GeneID','Editor']]
    df2_selected.rename(columns={'sgRNA sequence': 'Guide'}, inplace=True)
    df2_selected['Editor'] = df2_selected['Editor'].replace({'A-G': 'ABE', 'C-T': 'CBE'})
    table = merged_df.merge(df2_selected, on=["Guide","Editor"], how="left").drop_duplicates()
    
    new_order = ['Guide', 'GeneID', 'Editor'] + [col for col in table.columns if col not in ["Guide", "GeneID", "Editor"]]
    df_final = table[new_order]
    # df_final= df_final.sort_values(by="GeneID")
    # df_final = df_final.reset_index(drop=True)
    # df_final.rename(columns={'pDNA_pDNA_pDNA': 'pDNA'}, inplace=True)
    df_final.rename(columns={'GeneID': 'Gene'}, inplace=True)
    return df_final

#MC

def process_MC(df):
    df_melted = df.melt(id_vars=["guide", "Gene", "plasmid_reads"], var_name="Condition", value_name="Reads")

    df_melted["Cell_Line"] = df_melted["Condition"].str.extract(r'^(H23|HT29|PC9)')
    df_melted["Editor"] = df_melted["Condition"].str.extract(r'_(CBE|ABE)_')
    df_melted["Experiment"] = df_melted["Condition"].str.extract(r'_(T0|Sotor|Control|Adag|Tram|Pict|DebCet|Gefit|Osim)_')
    df_melted["Replicate"] = df_melted["Condition"].str.extract(r'_(exp\d+)_')

    df_melted = df_melted.drop(columns=["Condition"])

    df_melted["Condition"] = df_melted["Cell_Line"] + "_" + df_melted["Experiment"] + "_" + df_melted["Replicate"]

    df_melted = df_melted.drop(columns=["Cell_Line", "Experiment", "Replicate"])

    df_melted = df_melted.groupby(["guide", "Gene", "Editor", "plasmid_reads", "Condition"], as_index=False).sum()

    df_pivot = df_melted.pivot(index=["guide", "Gene", "Editor", "plasmid_reads"], 
                                columns="Condition", values="Reads").reset_index()

    df_pivot.columns.name = None

    df_pivot.columns = df_pivot.columns.str.replace(r'^HT29', 'HT_29', regex=True)
    df_pivot.columns = df_pivot.columns.str.replace(r'^PC9', 'PC_9', regex=True)
    df_pivot.columns = df_pivot.columns.str.replace(r'^H23', 'H_23', regex=True)
    #df_pivot = df_pivot.fillna("Na")
    df_pivot.rename(columns={'plasmid_reads': 'pDNA'}, inplace=True)
    df_pivot.rename(columns={'guide': 'Guide'}, inplace=True) #alterei isto
    df_pivot.columns = [re.sub(r'exp1', 'RepA', col) for col in df_pivot.columns]
    df_pivot.columns = [re.sub(r'exp2', 'RepB', col) for col in df_pivot.columns]
    df_pivot.columns = [re.sub(r'exp3', 'RepC', col) for col in df_pivot.columns]
    df_pivot.columns = df_pivot.columns.str.replace("__", "_", regex=False)

    MC_QC=pd.read_excel('data/Base Editing Screens Samplesheet.xlsx',sheet_name='BE sgRNA library_MC')
    MC_QC=MC_QC[['GeneID','sgRNA sequence','Editor']]
    MC_QC.rename(columns={'sgRNA sequence': 'Guide','GeneID':'Gene'}, inplace=True)

    MC_QC["Gene"] = MC_QC["Gene"].fillna("Unknown")
    MC_QC['Editor']=MC_QC['Editor'].replace({'A-G':'ABE','C-T':'CBE'})
    #only want the ones where the controls are
    MC_QC=MC_QC[MC_QC['Gene'].isin(['Panlethal splice donor', 'splice_nonessential',
        'Non-targeting control'])]

    MC=df_pivot.merge(MC_QC, on=["Guide","Editor"],how='left')
  

    MC["Gene_x"] = MC.apply(
            lambda row: row["Gene_y"] if pd.notna(row["Gene_y"]) else row["Gene_x"], 
            axis=1)

    MC = MC.drop(columns=["Gene_y"])
    MC.rename(columns={'Gene_x': 'Gene'}, inplace=True)
    return MC

def reorder_columns(df):
    cols = df.columns.tolist()
    new_order = ['Guide', 'Gene','Editor'] + [col for col in cols if col not in ['Guide', 'Gene','Editor']]
    return df[new_order]

def collapse_replicates_min(df):
    rep_dict = {}

    for col in df.columns:
        match = re.match(r'(.*)_(Rep[A-Z])$', col)
        if match:
            condition = match.group(1)
            rep = match.group(2)
            rep_dict.setdefault(condition, {})[rep] = col

    meta_cols = [col for col in ['Guide', 'Gene', 'Editor'] if col in df.columns]
    collapsed_df = df[meta_cols].copy()

    for condition, reps in rep_dict.items():
        rep_cols = list(reps.values())

        if len(rep_cols) >= 2:
            collapsed_df[condition] = df[rep_cols].apply(
                lambda row: row.iloc[row.abs().argmin()],
                axis=1
            )
        elif len(rep_cols) == 1:
            collapsed_df[condition] = df[rep_cols[0]]

    return collapsed_df

# Comparison names:
# RH_BRCA1_processed
# RH_BRCA2_processed
# RH_MCL1_processed
# RH_BCL2L1_processed
# RH_PARP1_processed
# MC_processed
# MC_processed_p
# GR_A549_ABE_CBE_MELJUSO_CBE_Pro
# GR_MELJUSO_ABE_Processed
# EG_pDNA
# EG_DO

def LFC_Z(table1,comparison_name):
    comparisons= pd.read_excel("data/Base Editing Screens Samplesheet.xlsx",sheet_name='Comparisons' )
    if comparison_name in comparisons['ComparisonName'].unique():
        comparisons= comparisons[comparisons['ComparisonName'] == comparison_name]    
    
    df_lfc = table1[['Guide', 'Gene', 'Editor']].copy() 


    for _, row in comparisons.iterrows():
        num = row['ConditionNumerator']
        denom = row['ConditionDenominator']

        if num in table1.columns and denom in table1.columns:
            table1[num] = pd.to_numeric(table1[num], errors='coerce')
            table1[denom] = pd.to_numeric(table1[denom], errors='coerce')           
            df_lfc[num] = np.log2((table1[num] +1) / (table1[denom]+1))  

    table1.iloc[:, 3:] = table1.iloc[:, 3:].apply(pd.to_numeric, errors='coerce')
    df_zscores=table1.copy()
    df_zscores = df_lfc.copy()
    
    
    lfc_columns = comparisons['ConditionNumerator'].unique()
    for col in lfc_columns:
        if col in df_lfc.columns:
            mean = df_lfc[col].mean()
            std = df_lfc[col].std()
            df_zscores[col] = (df_lfc[col] - mean) / std 
 
    return df_lfc, df_zscores


def process_EG(link1,link2,comparison):
    A375_DO=pd.read_csv(link1)
    A375_LIN=pd.read_csv(link2)
    A375_DO=A375_DO.drop(columns='Unnamed: 0')
    A375_LIN=A375_LIN.drop(columns='Unnamed: 0')
    HT29_DO=pd.read_csv('data/3_Counts/EG/counts_with_mismatches_CBE_activity_HT29_DO_PIC.csv')
    HT29_LIN=pd.read_csv('data/3_Counts/EG/counts_with_mismatches_CBE_activity_HT29_LIN_SCH.csv')
    HT29_DO=HT29_DO.drop(columns='Unnamed: 0')
    HT29_LIN=HT29_LIN.drop(columns='Unnamed: 0')
    A375_DO=A375_DO[['Construct Barcode','A375_RDA270_CP2165_PIC_RepA','A375_RDA270_CP2165_DO_RepB','A375_RDA270_CP2165_DO_RepA','A375_RDA270_CP2165_PIC_RepB','pDNA_CP2165']]
    A375_LIN=A375_LIN[['Construct Barcode','A375_RDA270_CP2165_SCH_RepB','pDNA_CP2165','A375_RDA270_CP2165_LIN_RepB','A375_RDA270_CP2165_SCH_RepA','A375_RDA270_CP2165_LIN_RepA']]
    A375_merge=A375_DO.merge(A375_LIN, on='Construct Barcode',suffixes=('_A375_DO_PIC','_A375_SCH_LIN'))
    A375_merge = A375_merge.rename(columns={'Construct Barcode': 'Guide'})
    A375_merge['Editor'] = 'CBE'
    HT29_DO=HT29_DO[['Construct Barcode','HT29_RDA270_CP2165_PIC_RepA','HT29_RDA270_CP2165_DO_RepB','HT29_RDA270_CP2165_DO_RepA','HT29_RDA270_CP2165_PIC_RepB','pDNA_CP2165']]
    HT29_LIN=HT29_LIN[['Construct Barcode','HT29_RDA270_CP2165_SCH_RepB','pDNA_CP2165','HT29_RDA270_CP2165_LIN_RepB','HT29_RDA270_CP2165_SCH_RepA','HT29_RDA270_CP2165_LIN_RepA']]
    HT29_merge=HT29_DO.merge(HT29_LIN, on='Construct Barcode',suffixes=('_HT29_DO_PIC','_HT29_SCH_LIN'))
    HT29_merge = HT29_merge.rename(columns={'Construct Barcode': 'Guide'})
    HT29_merge['Editor'] = 'CBE'
    assert A375_merge.shape[0]==A375_merge['Guide'].nunique()
    assert HT29_merge.shape[0]==HT29_merge['Guide'].nunique()
    EG_merge=HT29_merge.merge(A375_merge, on='Guide',suffixes=('_DO_PIC','_SCH_LIN'))
    EG_merge=EG_merge.drop(columns='Editor_DO_PIC')
    EG_merge=EG_merge.rename(columns={'Editor_SCH_LIN':'Editor'})
    processed=EG_merge.copy()
    bes_EG=pd.read_excel('data/Base Editing Screens Samplesheet.xlsx',sheet_name='BE sgRNA library_EG')
    bes_EG['Editor']=bes_EG['Editor'].replace({'A-G':'ABE','C-T':'CBE'})
    bes_EG=bes_EG[['sgRNA sequence', 'Editor', 'GeneID']].rename(columns={'sgRNA sequence':'Guide','GeneID':'Gene'})

    EG_merge=EG_merge.merge(bes_EG, on=['Guide','Editor'], how='left')
    EG_merge = reorder_columns(EG_merge)
    if comparison=='EG_DO':
        df_lfc, df_zscores=LFC_Z(EG_merge,comparison)
        df_lfc.to_csv("data/5_LFC/EG/EG_LFC_rep_pDNA.csv", index=False)
        df_zscores.to_csv("data/4_Screen_zscores/EG/EG_zscores_rep_pDNA.csv", index=False)
        zscores_min=collapse_replicates_min(df_zscores)
        zscores_min.to_csv("data/4_Screen_zscores/EG/EG_zscores_min_pDNA.csv", index=False)

    if comparison=='EG_pDNA': 
        df_lfc, df_zscores=LFC_Z(EG_merge,comparison)
        df_lfc.to_csv("data/5_LFC/EG/EG_LFC_rep_control.csv", index=False)
        df_zscores.to_csv("data/4_Screen_zscores/EG/EG_zscores_rep_control.csv", index=False)
        zscores_min=collapse_replicates_min(df_zscores)
        zscores_min.to_csv("data/4_Screen_zscores/EG/EG_zscores_min_control.csv", index=False)
    return processed,df_lfc,df_zscores,zscores_min






def process_df(file1,study,file2=None,sheet_name=None,comparison=None):
    zscores_min=None
    valid_studies = {'MC', 'GR', 'EG', 'RH'}
    valid_comparisons = {'RH_BRCA1_processed','RH_BRCA2_processed','RH_MCL1_processed',
        'RH_BCL2L1_processed','RH_PARP1_processed','MC_processed','MC_processed_p','GR_A549_ABE_CBE_MELJUSO_CBE_Pro',
        'GR_MELJUSO_ABE_Processed', 'EG_pDNA','EG_DO'}

    assert study in valid_studies, f"Invalid study: {study}. Must be one of {valid_studies}"
    if comparison is not None:
        assert comparison in valid_comparisons, f"Invalid comparison: {comparison}. Must be one of {valid_comparisons}"

    
    if study=='MC':
        df=pd.read_csv(f'data/3_Counts/MC/{file1}')
        df.sort_values(by='guide')
        df=process_MC(df)
        df.to_csv("data/3_Counts/MC/MC_processed.csv", index=False)
        if comparison=='MC_processed_p':
            lfc_data,zscores=LFC_Z(df,comparison)
            lfc_data.to_csv("data/5_LFC/MC/MC_LFC_rep_pDNA.csv", index=False)
            zscores.to_csv("data/4_Screen_zscores/MC/MC_zscores_rep_pDNA.csv", index=False)
            zscores_min=collapse_replicates_min(zscores)
            zscores_min.to_csv("data/4_Screen_zscores/MC/MC_zscores_min_pDNA.csv", index=False)
        if comparison=='MC_processed':
            lfc_data,zscores=LFC_Z(df,comparison)
            lfc_data.to_csv("data/5_LFC/MC/MC_LFC_rep_control.csv", index=False)
            zscores.to_csv("data/4_Screen_zscores/MC/MC_zscores_rep_control.csv", index=False)
            zscores_min=collapse_replicates_min(zscores)
            zscores_min.to_csv("data/4_Screen_zscores/MC/MC_zscores_min_control.csv", index=False)
    if study=='EG':
        df,lfc_data,zscores,zscores_min=process_EG(file1,file2,comparison)

        
    if study=='GR':
        name=file1.replace(r'.csv', '').replace(r'.tsv', '').replace(r'.txt', '')
        df=pd.read_csv(f"data/3_Counts/GR/{file1}")
        df.sort_values(by="Construct Barcode") #all unique guides
        df=process_GR(df)
        df.to_csv(f"data/3_Counts/GR/GR_{name}_processed.csv", index=False)
        if comparison== 'GR_A549_ABE_CBE_MELJUSO_CBE_Pro':
            lfc_data, zscores=LFC_Z(df,"GR_A549_ABE_CBE_MELJUSO_CBE_Pro")
            lfc_data.to_csv("data/5_lfc/GR/GR_A549_ABE_CBE_MELJUSO_CBE_LFC.csv", index=False)
            zscores.to_csv("data/4_Screen_zscores/GR/GR_A549_ABE_CBE_MELJUSO_CBE_zscores.csv", index=False)
        if comparison=='GR_MELJUSO_ABE_Processed': 
            lfc_data, zscores=LFC_Z(df,"GR_MELJUSO_ABE_Processed")
            lfc_data.to_csv("data/5_lfc/GR/GR_MELJUSO_ABE_LFC.csv", index=False)
            zscores.to_csv("data/4_Screen_zscores/GR/GR_MELJUSO_ABE_zscores.csv", index=False)

    if study=='RH':
        if sheet_name=='BE3_HAP1':
            subset='PARP1_set'
            df=pd.read_excel(f"data/3_Counts/RH/{file1}",header=None, sheet_name=sheet_name)
            df=process_RH(df,subset=subset)
            df.to_csv("data/3_Counts/RH/RH_PARP1_processed.csv", index=False)
            lfc_data,zscores=LFC_Z(df,'RH_PARP1_processed')
            lfc_data.to_csv("data/5_LFC/RH/RH_PARP1_LFC.csv", index=False)
            zscores.to_csv("data/4_Screen_zscores/RH/RH_PARP1_zscores.csv", index=False)

        if sheet_name=='BRCA1' or sheet_name=='BRCA2':
            subset='BRCA1_BRCA2_set'
            df=pd.read_excel(f"data/3_Counts/RH/{file1}",header=None, sheet_name=sheet_name)
            df=process_RH(df,subset=subset)
            df.to_csv(f"data/3_Counts/RH/RH_{sheet_name}_processed.csv", index=False)
            if sheet_name=='BRCA1':
                lfc_data, zscores=LFC_Z(df,"RH_BRCA1_processed")
                lfc_data.to_csv("data/5_LFC/RH/RH_BRCA1_LFC.csv", index=False)
                zscores.to_csv("data/4_Screen_zscores/RH/RH_BRCA1_zscores.csv", index=False)
            if sheet_name=='BRCA2':
                lfc_data, zscores=LFC_Z(df,"RH_BRCA2_processed")
                lfc_data.to_csv("data/5_LFC/RH/RH_BRCA2_LFC.csv", index=False)
                zscores.to_csv("data/4_Screen_zscores/RH/RH_BRCA2_zscores.csv", index=False)
                                
        if sheet_name=='MCL1' or sheet_name=='BCL2L1':
            subset='BCL2L1_MCL1_set'
            df=pd.read_excel(f"data/3_Counts/RH/{file1}",header=None, sheet_name=sheet_name)
            df=process_RH(df,subset=subset)
            df.to_csv(f"data/3_Counts/RH/RH_{sheet_name}_processed.csv", index=False)
            lfc_data, zscores=LFC_Z(df,"RH_MCL1_processed")
            lfc_data.to_csv("data/5_LFC/RH/RH_MCL1_LFC.csv", index=False)
            zscores.to_csv("data/4_Screen_zscores/RH/RH_MCL1_zscores.csv", index=False)
    return df,lfc_data,zscores , zscores_min
    

    


def process_all(df,study,subset):
    if study=="RH":
        processed=process_RH(df,subset=subset)
        df_lfc, df_zscores = LFC_Z(processed)
    elif study=="GR":
        processed=process_GR(df)
        df_lfc, df_zscores = LFC_Z(processed)
    elif study=="MC":
        processed=process_MC(df)
        df_lfc, df_zscores = LFC_Z(processed)
    else:
        print('No processing for that study yet')
    return processed, df_lfc, df_zscores

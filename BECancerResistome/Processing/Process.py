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

    table=df.merge(df2_selected, on="Guide",how="left").drop_duplicates()
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
    df2_selected = bes[['sgRNA sequence','GeneID']]
    df2_selected.rename(columns={'sgRNA sequence': 'Guide'}, inplace=True)
    # df2_selected['Editor'] = df2_selected['Editor'].replace({'A-G': 'ABE', 'C-T': 'CBE'})
    table = merged_df.merge(df2_selected, on=["Guide"], how="left").drop_duplicates()
    #table=merged_df.merge(df2_selected, on=["Guide","Editor"],how="left").drop_duplicates()
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

    df_pivot.columns = df_pivot.columns.astype(str).str.replace(r"(?<!^)(?=[A-Z])", "_", regex=True)
    df_pivot = df_pivot.fillna("Na")
    df_pivot.rename(columns={'plasmid_reads': 'pDNA'}, inplace=True)
    df_pivot.rename(columns={'guide': 'Guide'}, inplace=True) #alterei isto
    df_pivot.columns = [re.sub(r'exp1', 'RepA', col) for col in df_pivot.columns]
    df_pivot.columns = [re.sub(r'exp2', 'RepB', col) for col in df_pivot.columns]
    df_pivot.columns = [re.sub(r'exp3', 'RepC', col) for col in df_pivot.columns]
    df_pivot.columns = df_pivot.columns.str.replace("__", "_", regex=False)

    MC_QC=pd.read_excel('data/41588_2024_1948_MOESM4_ESM.xlsx',sheet_name='ST1 base_editing')
    MC_QC=MC_QC[['Gene','guide','sgRNA_type']]
    MC_QC.rename(columns={'guide': 'Guide'}, inplace=True)
    MC_QC["Gene"] = MC_QC["Gene"].fillna("Unknown")
    print(MC_QC['sgRNA_type'].unique())


    MC=df_pivot.merge(MC_QC, on=["Guide","Gene"],how='left')
    #Replace controls in Gene column
    MC["Gene"] = MC.apply(lambda row: row["sgRNA_type"] if row["sgRNA_type"] != "exonic" else row["Gene"], axis=1)
    MC = MC.drop(columns=["sgRNA_type"])
    return(MC)


# Comparison names:
# RH_BRCA1_processed
# RH_BRCA2_processed
# RH_MCL1_processed
# RH_BCL2L1_processed
# RH_PARP1_processed
# MC_processed
# GR_A549_ABE_CBE_MELJUSO_CBE_Pro
# GR_MELJUSO_ABE_Processed

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

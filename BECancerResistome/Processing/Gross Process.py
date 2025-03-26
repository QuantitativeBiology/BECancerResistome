import pandas as pd
import numpy as np
import re

#RH

def process_RH(df,subset='BRCA1_BRCA2_set'):

    # Combine the first three rows into one header
    new_header = df.iloc[:3].apply(lambda x: '_'.join(x), axis=0)
    df.columns = new_header
    df = df.drop(index=[0, 1, 2]) # Drop the first three rows since they're now part of the header

    df.reset_index(drop=True, inplace=True)
    df.rename(columns={'Cell line_Condition_sgRNA sequence': 'Guide'}, inplace=True)

    bes= pd.read_excel("/Users/joanacorreia/Desktop/Tese/Original Counts files/Base Editing Screens Samplesheet.xlsx",sheet_name='BE sgRNA library_RH' )
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

def RH_LFC_Z(table1):
    reference_col = "pDNA"

    table1.iloc[:, 3:] = table1.iloc[:, 3:].apply(pd.to_numeric) #just so no errors
    table1[reference_col] = table1[reference_col].astype(float)
    
    df_lfc = table1.copy()
    df_zscores=table1.copy()
    for i in range(3,(table1.shape)[1]):
        df_lfc.iloc[:, i]= np.log2(pd.to_numeric((table1.iloc[:, i] + 1) / (table1[reference_col][:] + 1)))
        if i!=1:
            df_zscores.iloc[:, i]= (df_lfc.iloc[:, i] - np.mean(df_lfc.iloc[:, i])) / np.std(df_lfc.iloc[:, i])
        else:
            df_zscores.iloc[:, i]= 0
    # df_zscores =df_zscores.drop["pDNA"]
    # df_lfc= df_lfc.drop["pDNA"]
    return df_lfc, df_zscores


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
    bes= pd.read_excel("/Users/joanacorreia/Desktop/Tese/Original Counts files/Base Editing Screens Samplesheet.xlsx",sheet_name='BE sgRNA library_GR' )
    df2_selected = bes[['sgRNA sequence', 'Editor','GeneID']]
    df2_selected.rename(columns={'sgRNA sequence': 'Guide'}, inplace=True)
    df2_selected['Editor'] = df2_selected['Editor'].replace({'A-G': 'ABE', 'C-T': 'CBE'})

    table=merged_df.merge(df2_selected, on=["Guide","Editor"],how="left").drop_duplicates()
    new_order = ['Guide', 'GeneID', 'Editor'] + [col for col in table.columns if col not in ["Guide", "GeneID", "Editor"]]
    df_final = table[new_order]
    df_final= df_final.sort_values(by="GeneID")
    df_final = df_final.reset_index(drop=True)
    df_final.rename(columns={'pDNA_pDNA_pDNA': 'pDNA'}, inplace=True)
    df_final.rename(columns={'GeneID': 'Gene'}, inplace=True)

    return df_final


def GR_LFC_Z(table1):
    reference_col = "pDNA"

    table1.iloc[:, 3:] = table1.iloc[:, 3:].apply(pd.to_numeric) #just so no errors

    table1[reference_col] = table1[reference_col].astype(float)
    df_lfc = table1.copy()
    df_zscores=table1.copy()
    for i in range(3,(table1.shape)[1]):
        df_lfc.iloc[:, i]= np.log2(pd.to_numeric((table1.iloc[:, i] + 1) / (table1[reference_col][:] + 1)))
        if i!=1:
            df_zscores.iloc[:, i]= (df_lfc.iloc[:, i] - np.mean(df_lfc.iloc[:, i])) / np.std(df_lfc.iloc[:, i])
        else:
            df_zscores.iloc[:, i]= 0
    # df_zscores =df_zscores.drop["pDNA"]
    # df_lfc= df_lfc.drop["pDNA"]
    return df_lfc, df_zscores



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
    df_pivot.columns = [re.sub(r'exp1', 'repA', col) for col in df_pivot.columns]
    df_pivot.columns = [re.sub(r'exp2', 'repB', col) for col in df_pivot.columns]
    df_pivot.columns = [re.sub(r'exp3', 'repC', col) for col in df_pivot.columns]
    return(df_pivot)

def MC_LFC_Z(a):
    a.columns = [re.sub(r'repA', 'exp1',col) for col in a.columns]
    a.columns = [re.sub(r'repB','exp2', col) for col in a.columns]
    a.columns = [re.sub(r'repC','exp3', col) for col in a.columns]

    control_cols = {col: re.search(r'exp\d+', col).group() for col in a.columns if '__Control_exp' in col}

    non_experimental_cols = ['Guide', 'Gene','Editor']
    lfc_data = a[non_experimental_cols].copy() 

    # Compute LFC for each condition
    for col in a.columns:
        if '__Control_exp' not in col and col not in non_experimental_cols:  #
            match = re.search(r'exp\d+', col)  # Extract experiment number 1 or 2
            if match:
                exp_num = match.group()  #like, "exp1"
                
                # Find the corresponding control column
                control_col = next((ctrl for ctrl, exp in control_cols.items() if exp == exp_num), None)
                
                if control_col:
                    lfc_data[col] = np.log2((a[col]+1 )/ (a[control_col]+1))  

    lfc_data.columns = [re.sub(r'exp1', 'repA', col) for col in lfc_data.columns]
    lfc_data.columns = [re.sub(r'exp2', 'repB', col) for col in lfc_data.columns]
    lfc_data.columns = [re.sub(r'exp3', 'repC', col) for col in lfc_data.columns]


   

    df_zscores = lfc_data.copy()  

    for col in df_zscores.columns:
        if col not in non_experimental_cols:  
            mean = df_zscores[col].mean()
            std = df_zscores[col].std()
            
            if std == 0:  
                df_zscores[col] = 0
            else:
                df_zscores[col] = (df_zscores[col] - mean) / std

    return lfc_data,df_zscores



def process_all(df,study,subset):
    if study=="RH":
        processed=process_RH(df,subset=subset)
        df_lfc, df_zscores = RH_LFC_Z(processed)
    elif study=="GR":
        processed=process_GR(df)
        df_lfc, df_zscores = GR_LFC_Z(processed)
    elif study=="MC":
        processed=process_MC(df)
        df_lfc, df_zscores = MC_LFC_Z(processed)
    else:
        print('No processing for that study yet')
    return processed, df_lfc, df_zscores

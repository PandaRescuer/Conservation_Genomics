import time
import pandas as pd

f_LOF = pd.read_csv('LOF_snp_RMA_SFS_no_Miss.txt', sep='\t')
f_Del = pd.read_csv('mis_snp_DELETERIOUS_SFS_no_Miss.txt', sep='\t')
f_Tol = pd.read_csv('mis_snp_TOLERATED_SFS_no_Miss.txt', sep='\t')
f_Mis = pd.read_csv('missense_snp_RMA_SFS_no_Miss.txt', sep='\t')
f_Syn = pd.read_csv('synonymous_snp_RMA_SFS_no_Miss.txt', sep='\t')
f_CNE = pd.read_csv('CNE_snp_RMA_SFS_no_Miss.txt', sep='\t')

# Merge all the data
merged_data = pd.concat([f_CNE['num'], f_Del['num'], f_LOF['num'],f_Mis['num'],  f_Syn['num'], f_Tol['num']], axis=1)

# Rename the merged columns
merged_data.columns = ['CNE', 'Deleterous', 'LOF','Missense', 'Synonymous', 'Tolerated']

# Add columns a and b from the original data to the merged data
merged_data['POP'] = f_LOF['POP']
merged_data['freq'] = f_LOF['freq']
merged_data["DAF"] = (
    merged_data["freq"]
    / merged_data.groupby("POP")["freq"].transform("max")
)

# Reorder the columns if needed
merged_data = merged_data[['POP', 'freq', 'DAF', 'CNE', 'Deleterous', 'LOF','Missense', 'Synonymous', 'Tolerated']]

# Output the merged result
print(merged_data)
merged_data.to_csv('merge_SFS_no_Miss.txt', sep='\t', index=False)

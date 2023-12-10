import pandas as pd

# Define column names
column_names = ['SW_score', 'perc_divergence', 'perc_deletion', 'perc_insertion', 'query_sequence', 
                'query_begin', 'query_end', 'query_remaining', 'strand', 'matching_repeat', 
                'repeat_class/family', 'repeat_begin', 'repeat_end', 'repeat_remaining', 'ID', 'star']

# Read the first few lines of your data
df = pd.read_csv("recentintrons_multidnaseq.fa.out", sep='\s+', names=column_names, skiprows=3, header=None)

print(df.head())


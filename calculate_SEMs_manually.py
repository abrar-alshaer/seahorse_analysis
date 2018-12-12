#calculating SEMs manually from pandas dataframe using groupby to get groups of biological replicates then using groupby.size() to get n (sample size)
''' Calculating SEMs manually
print "\nSIZE (n) OF EACH GROUP:", df.groupby('2) Group').size(), "\n"
sem_bio_reps = []
bioGroups = []
count = 0
for index, row in std_bio_reps.iterrows():
    bioGroups.append(row[0])
    sem_bio_reps.append(row[1]/math.sqrt(df.groupby('2) Group').size()[count]))
    print row[0], row[1]/math.sqrt(df.groupby('2) Group').size()[count])
    count = count+1
'''

#helpful code for transposing and maintaining row/column names in place
#transposing all the dataframes 
# avg_t = avg_bio_reps.set_index('2) Group').T
# std_t = std_bio_reps.set_index('2) Group').T
# sem_t = sem_bio_reps.set_index('2) Group').T

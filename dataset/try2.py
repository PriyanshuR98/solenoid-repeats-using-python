import scipy.stats as stats
import matplotlib.pyplot as plt
from collections import Counter




dicfact1={}
dicfact1={'A':-0.591 , 'C':-1.343,'D':1.050 , 'E':1.375 , 'F':-1.006 , 'G':-0.384 , 'H':0.036 ,'I':-1.29,'K':1.831,'L':-1.019,'M':-0.663,'N':0.945,'P':0.189,'Q':0.931,'R':1.538,'S':-0.228,'T':-0.032,'V':-1.337,'W':-0.595,'Y':0.260}
s="GAMGISLGNSEADRQLLEAAKAGDVETVKKLCTVQSVNCRDIEGRQSTPLHFAAGYNRVSVVEYLLQHGADVHAKDKGGLVPLHNACSYGHYEVAELLVKHGAVVNVADLWKFTPLHEAAAKGKYEICKLLLQHGADPTKKNRDGNTPLDLVKDGDTDIQDLLRGDAALLDAAKK"
lenzz=(len(s))
ans1=""
vec1=list();
for i in range(len(s)):
    ans1+=str(dicfact1[s[i]])+" "
    vec1.append((dicfact1[s[i]]))

values = []
j=0
for i in vec1:
    values.append(i);




zscores = stats.zscore(values)
print(zscores)


# frequency_rank = []
# for i in range(len(s)):
#     amino_acids_at_position = [sequence[i] for sequence in s]
#     amino_acid_counts = Counter(amino_acids_at_position)
#     rank = {amino_acid: count/len(s) for amino_acid, count in amino_acid_counts.items()}
#     frequency_rank.append(rank)



plt.plot(zscores)
plt.show()

# # print(*zz)
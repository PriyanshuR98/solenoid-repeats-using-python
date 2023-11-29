import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
from collections import Counter
import scipy.stats as stats


def calculate_stddev(numbers):
    mean = sum(numbers) / len(numbers)
    sum_of_squared_diff = sum((x.real - mean) ** 2 for x in numbers)
    variance = sum_of_squared_diff / len(numbers)
    stddev = abs(variance)
    return stddev


def calc_avg(num):
    mean = sum(num)/len(num)
    return mean


def z_score(protein_seq, dft_mean, dft_std):
    dft = np.sum(
        [1 if aa in 'DE' else -1 if aa in 'RK' else 0 for aa in protein_seq])
    z = (dft - dft_mean) / dft_std
    return z


def zScore(vec1):
    zscV1 = list()
    average = calc_avg(vec1)
    stdDeviation = calculate_stddev(vec1)
    for i in range(0, len(vec1)):
        zscV1[i] = (vec1[i]-average)/(stdDeviation)
    return zscV1


def calculate_mode(numbers):
    frequency = Counter(numbers)
    max_frequency = max(frequency.values())
    modes = [num for num, freq in frequency.items() if freq == max_frequency]
    return modes


def compute_theta_ratio(metric_functions, Nsp):
    theta_ratios = []
    for theta in range(1, 6):
        N_theta = 0

        # Loop over the five metric functions (a = 1, ..., 5)
        for a in range(0, len(metric_functions)):
            for n in range(1, Nsp):
                zn_a = metric_functions[a][n]
                if zn_a > theta:
                    N_theta += 1

        # print(N_theta)
        theta_ratio = (100 * N_theta) / (5 * (Nsp - 1))
        theta_ratios.append(theta_ratio)
        

    return theta_ratios


dicfact1 = {}
dicfact1 = {'A': -0.591, 'C': -1.343, 'D': 1.050, 'E': 1.375, 'F': -1.006, 'G': -0.384, 'H': 0.036, 'I': -1.29, 'K': 1.831, 'L': -
            1.019, 'M': -0.663, 'N': 0.945, 'P': 0.189, 'Q': 0.931, 'R': 1.538, 'S': -0.228, 'T': -0.032, 'V': -1.337, 'W': -0.595, 'Y': 0.260}

dicfact2 = {}
dicfact2 = {'A': -1.302, 'C': 0.465, 'D': 0.302, 'E': -1.453, 'F': -0.590, 'G': 1.652, 'H': -0.417, 'I': -0.547, 'K': -0.561, 'L': -
            0.987, 'M': -1.524, 'N': 0.828, 'P': 2.081, 'Q': -0.179, 'R': -0.055, 'S': 1.399, 'T': 0.326, 'V': -0.279, 'W': 0.009, 'Y': 0.830}

dicfact3 = {}
dicfact3 = {'A': -0.733, 'C': -0.862, 'D': -3.656, 'E': 1.477, 'F': 1.891, 'G': 1.330, 'H': -1.673, 'I': 2.131, 'K': 0.533, 'L': -
            1.505, 'M': 2.219, 'N': 1.299, 'P': -1.628, 'Q': -3.005, 'R': 1.502, 'S': -4.760, 'T': 2.213, 'V': -0.544, 'W': 0.672, 'Y': 3.097}

dicfact4 = {}
dicfact4 = {'A': 1.570, 'C': -1.020, 'D': -0.259, 'E': 0.113, 'F': -0.397, 'G': 1.045, 'H': -1.474, 'I': 0.393, 'K': -0.277, 'L': 1.266,
            'M': -1.005, 'N': -0.169, 'P': 0.421, 'Q': -0.503, 'R': 0.440, 'S': 0.670, 'T': 0.908, 'V': 1.242, 'W': -2.128, 'Y': -0.838}

dicfact5 = {}
dicfact5 = {'A': -0.146, 'C': -0.255, 'D': -3.242, 'E': -0.837, 'F': 0.412, 'G': 2.064, 'H': -0.078, 'I': 0.816, 'K': 1.648, 'L': -
            0.912, 'M': 1.212, 'N': 0.933, 'P': -1.392, 'Q': -1.853, 'R': 2.897, 'S': -2.647, 'T': 1.313, 'V': -1.262, 'W': -0.184, 'Y': 1.512}


# The repeating protein sequences that here is taken (3TWQ|Chains A, B|Tankyrase-2|Homo sapiens (9606)) from PDB
# s="GAMGISLGNSEADRQLLEAAKAGDVETVKKLCTVQSVNCRDIEGRQSTPLHFAAGYNRVSVVEYLLQHGADVHAKDKGGLVPLHNACSYGHYEVAELLVKHGAVVNVADLWKFTPLHEAAAKGKYEICKLLLQHGADPTKKNRDGNTPLDLVKDGDTDIQDLLRGDAALLDAAKK"

# The Repeating sequnece mentioned in reasearch paper of ID:-1EZG
s = "QCTGGADCTSCTGACTGCGNCPNAVTCTNSQHCVKANTCTGSTDCNTAQTCTNSKDCFEANTCTDSTNCYKATACTNSSGCPGH"

# print(len(s))
ans1 = ""
vec1 = list()
for i in range(len(s)):
    ans1 += str(dicfact1[s[i]])+" "
    vec1.append((dicfact1[s[i]]))

ans2 = ""
vec2 = list()
for i in range(len(s)):
    ans2 += str(dicfact2[s[i]])+" "
    vec2.append((dicfact2[s[i]]))

ans3 = ""
vec3 = list()
for i in range(len(s)):
    ans3 += str(dicfact3[s[i]])+" "
    vec3.append((dicfact3[s[i]]))

ans4 = ""
vec4 = list()
for i in range(len(s)):
    ans4 += str(dicfact4[s[i]])+" "
    vec4.append((dicfact4[s[i]]))

ans5 = ""
vec5 = list()
for i in range(len(s)):
    ans5 += str(dicfact5[s[i]])+" "
    vec5.append((dicfact5[s[i]]))


vec1 = np.fft.fft(vec1)
vec2 = np.fft.fft(vec2)
vec3 = np.fft.fft(vec3)
vec4 = np.fft.fft(vec4)
vec5 = np.fft.fft(vec5)


lenx = len(s)

for i in range(1, lenx):

    a = np.real(vec1[i])
    b = np.imag(vec1[i])
    result = np.fabs((a**2) + (b**2))
    result = result / (math.sqrt(lenx))
    result = result/2
    # a=a/(np.fabs(lenx))
    vec1[i-1] = result

for i in range(1, lenx):
    a = np.real(vec2[i])
    b = np.imag(vec2[i])
    result = np.fabs((a**2) + (b**2))
    result = result / (math.sqrt(lenx))
    # a=a/(np.fabs(lenx))
    vec2[i-1] = result

for i in range(1, lenx):
    a = np.real(vec3[i])
    b = np.imag(vec3[i])
    result = np.fabs((a**2) + (b**2))
    result = result / (math.sqrt(lenx))
    # a=a/(np.fabs(lenx))
    vec3[i-1] = result

for i in range(1, lenx):
    a = np.real(vec4[i])
    b = np.imag(vec4[i])
    result = np.fabs((a**2) + (b**2))
    result = result / (math.sqrt(lenx))
    # a=a/(np.fabs(lenx))
    vec4[i-1] = result

for i in range(1, lenx):
    a = np.real(vec5[i])
    b = np.imag(vec5[i])
    result = np.fabs((a**2) + (b**2))
    result = result / (np.fabs(lenx))
    # a=a/(np.fabs(lenx))
    vec5[i-1] = result


zscores_vec1 = stats.zscore(vec1)
zscores_vec2 = stats.zscore(vec2)
zscores_vec3 = stats.zscore(vec3)
zscores_vec4 = stats.zscore(vec4)
zscores_vec5 = stats.zscore(vec5)

# zscores_vec1=zscores_vec1.view(dtype=np.real)
for i in range(0, lenx):
    zscores_vec1[i] = np.real(zscores_vec1[i])
for i in range(0, lenx):
    zscores_vec2[i] = np.real(zscores_vec2[i])
for i in range(0, lenx):
    zscores_vec3[i] = np.real(zscores_vec3[i])
for i in range(0, lenx):
    zscores_vec4[i] = np.real(zscores_vec4[i])
for i in range(0, lenx):
    zscores_vec5[i] = np.real(zscores_vec5[i])
# print(zscores_vec1)

maxiZscore = -1e9
maxiZscore = max(maxiZscore, max(zscores_vec1))
maxiZscore = max(maxiZscore, max(zscores_vec2))
maxiZscore = max(maxiZscore, max(zscores_vec3))
maxiZscore = max(maxiZscore, max(zscores_vec4))
maxiZscore = max(maxiZscore, max(zscores_vec5))

theta_ratios = compute_theta_ratio(
    [zscores_vec1, zscores_vec2, zscores_vec3, zscores_vec4, zscores_vec5], lenx)
# maxiZscore = max(zscores_vec1,zscores_vec2,zscores_vec3,zscores_vec4,zscores_vec5)
print('thehta rations', theta_ratios)
print('z score max', maxiZscore)

# Generate some random data for plotting
x = np.linspace(0, 10, 100)

# Create a figure and subplots
fig, axes = plt.subplots(2, 3, figsize=(12, 8))

# Plot the graphs on each subplot
axes[0, 0].plot(zscores_vec1)
axes[0, 0].set_title('Property 1')

axes[0, 1].plot(zscores_vec2)
axes[0, 1].set_title('Property 2')

axes[0, 2].plot(zscores_vec3)
axes[0, 2].set_title('Property 3')

axes[1, 0].plot(zscores_vec4)
axes[1, 0].set_title('Property 4')

axes[1, 1].plot(zscores_vec5)
axes[1, 1].set_title('Property 5')

# Hide the empty subplot
axes[1, 2].axis('off')

# Adjust the spacing between subplots
plt.tight_layout()

# Show the plot
plt.show()


# x1=np.array(vec1)
# op1=np.fft.fft(x1)
# len_op1=len(s)-1

# # for i in range(len(op1)):
# #     # op1[i]=(abs(op1[i]))
# #     op1[i]=op1[i].real
# # print(*op1)

# std_dev1=calculate_stddev(op1)
# average1=calc_avg(op1)
# average1=average1.real
# # print(average1)

# zz=[]
# for i in range(len(s)):
#     z_sore=z_score(s[i],average1,std_dev1)
#     zz.append(z_sore)

# res1=calculate_mode(zz)
# print(*zz)

# plot1=[]
# z_index=[]

# for i in range(len(zz)):
#     if(zz[i] >= res1):
#         plot1.append(op1[i])
#         z_index.append(i+1)


# N1=len(zz)
# plt.plot(z_index,plot1)
# plt.show()
# freq1=np.arange(N1)
# # plt.plot(freq, np.abs(op))
# # plt.show()
# x2=np.array(vec2)
# op2=np.fft.fft(x2)
# N2=len(op2)
# freq2=np.arange(N2)
# #plt.plot(freq, np.abs(op))
# #plt.show()
# x3=np.array(vec3)
# op3=np.fft.fft(x3)
# N3=len(op3)
# freq3=np.arange(N3)
# #plt.plot(freq, np.abs(op))
# #plt.show()
# x4=np.array(vec4)
# op4=np.fft.fft(x4)
# N4=len(op4)
# freq4=np.arange(N4)
# #plt.plot(freq, np.abs(op))
# #plt.show()
# x5=np.array(vec5)
# op5=np.fft.fft(x5)
# N5=len(op5)
# freq5=np.arange(N5)
# #plt.plot(freq, np.abs(op))
# #plt.show()

# fig,axis_plot=plt.subplots(nrows=5, ncols=1, figsize=(10, 12))
# axis_plot[0].plot(freq1, np.abs(op1))
# axis_plot[0].set_title('Factor 1', fontsize=8)
# axis_plot[1].plot(freq2, np.abs(op2))
# axis_plot[1].set_title('Factor 2', fontsize=8)
# axis_plot[2].plot(freq3, np.abs(op3))
# axis_plot[2].set_title('Factor 3', fontsize=8)
# axis_plot[3].plot(freq4, np.abs(op4))
# axis_plot[3].set_title('Factor 4', fontsize=8)
# axis_plot[4].plot(freq5, np.abs(op5))
# axis_plot[4].set_title('Factor 5', fontsize=8)


# # fig.suptitle('Protein Sequence DFT')
# fig.suptitle('Protein Sequence DFT', fontsize=8, y=0.02)
# # plt.tight_layout()
# plt.show()
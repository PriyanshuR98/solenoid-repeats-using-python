import pandas as pd
from Bio import *
from Bio.PDB import PDBList, PDBParser, PPBuilder
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
from collections import Counter
import scipy.stats as stats
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, matthews_corrcoef
import pickle
import requests
"""Required Functions"""

def compute_theta_ratio(metric_functions, Nsp):
    theta_ratios = []
    for theta in range(1, 6):
        N_theta = 0
        for a in range(0, len(metric_functions)):
            for n in range(1, Nsp):
                zn_a = metric_functions[a][n]
                if zn_a > theta:
                    N_theta += 1
        theta_ratio = (100 * N_theta) / (5 * (Nsp - 1))
        theta_ratios.append(theta_ratio)
    return theta_ratios


data = pd.read_csv('DATA.csv')
# column_0_values = data.iloc[0:200, 0].tolist()
# test_0_values = data.iloc[570:610, 0].tolist()
column_0_values = data.iloc[0:20, 0].tolist()
test_0_values = data.iloc[20:22, 0].tolist()
print(test_0_values)

# def get_protein_fasta_sequence(pdb_id):
#     pdb_list = PDBList()
#     temp_dir = "."

#     pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir=temp_dir,file_format="pdb")

#     pdb_parser = PDBParser(QUIET=True)

#     structure = pdb_parser.get_structure(pdb_id, pdb_filename)
#     pp_builder = PPBuilder()

#     sequences = []
#     for model in structure:
#         for chain in model:
#             peptides = pp_builder.build_peptides(chain)
#             for peptide in peptides:
#                 sequence = str(peptide.get_sequence())
#                 sequences.append(sequence)
#     import os
#     os.remove(pdb_filename)

#     return sequences

# # pdb_id = "1ezg"
# lst=[]
# for x in column_0_values:
#   protein_sequences = get_protein_fasta_sequence(x)
#   # print(protein_sequences)
#   lst.append(protein_sequences)

# # print(lst)
# seqMatrix=[]
# for x in lst:
#   for i in x:
#     seqMatrix.append(i)
#     break

# # print(matrix)
# # print(len(seqMatrix))

"""Conversion of PDBID to UNIPROTID"""


def fetch_uniprot_id_from_pdb(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        try:
            # Extract UniProt ID from the JSON response
            for entity in data['rcsb_entry_container_identifiers']['polymer_entity_ids']:
                entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity}"
                entity_response = requests.get(entity_url)
                if entity_response.status_code == 200:
                    entity_data = entity_response.json()
                    for reference in entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids', []):
                        return reference
            return "notFound"
        except KeyError:
            return "Unable to extract UniProt ID"
    else:
        return f"Error fetching data: HTTP {response.status_code}"



"""Fetching Repeats and Non Repeating Unit Information from UNIPROT API using UNIPROTID

"""

def fetch_uniprot_data(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)

    if response.status_code == 200:
        return response.text
    else:
        return f"Error fetching data: HTTP {response.status_code}"

# Function to Extract Repeat Information


def extract_repeat_info(uniprot_data):
    lines = uniprot_data.split("\n")
    repeat_info = [line for line in lines if line.startswith("FT   REPEAT")]
    return True if repeat_info else False


"""Presence of Repeats and Non Repeats"""

pdb = PDBList()
resultRepeatNR = []
resultantList = []
testResultantList = []
testList = []
print(column_0_values)
for x in column_0_values:
    uniprotID = fetch_uniprot_id_from_pdb(x)
    if uniprotID != "notFound":
        resultantList.append(x)
        checkFetchUNI = fetch_uniprot_data(uniprotID)
        checkRepeatFinalUNI = extract_repeat_info(checkFetchUNI)
        if checkRepeatFinalUNI:
            resultRepeatNR.append(1)
        else:
            resultRepeatNR.append(0)
    else:
        continue
    # check =
# print(resultRepeatNR)
# print(len(resultRepeatNR))

# for x in test_0_values:
#   uniprotID = fetch_uniprot_id_from_pdb(x)
#   if uniprotID != "notFound":
#     testResultantList.append(x)
#     checkFetchUNI = fetch_uniprot_data(uniprotID)
#     checkRepeatFinalUNI = extract_repeat_info(checkFetchUNI)
#     if checkRepeatFinalUNI:
#       testList.append(1)
#     else:
#       testList.append(0)
#   else:
#     continue

"""Testing Check"""

# for x in test_0_values:
#   uniprotID = fetch_uniprot_id_from_pdb(x)
#   if uniprotID != "notFound":
#     testResultantList.append(x)
#     checkFetchUNI = fetch_uniprot_data(uniprotID)
#     checkRepeatFinalUNI = extract_repeat_info(checkFetchUNI)
#     if checkRepeatFinalUNI:
#       testList.append(1)
#     else:
#       testList.append(0)
#   else:
#     continue

"""Getting sequences for the Required ones

"""


def get_protein_fasta_sequence(pdb_id):
    pdb_list = PDBList()
    temp_dir = "."

    pdb_filename = pdb_list.retrieve_pdb_file(
        pdb_id, pdir=temp_dir, file_format="pdb")

    pdb_parser = PDBParser(QUIET=True)

    structure = pdb_parser.get_structure(pdb_id, pdb_filename)
    pp_builder = PPBuilder()

    sequences = []
    for model in structure:
        for chain in model:
            peptides = pp_builder.build_peptides(chain)
            for peptide in peptides:
                sequence = str(peptide.get_sequence())
                sequences.append(sequence)
    import os
    os.remove(pdb_filename)

    return sequences


# pdb_id = "1ezg"
lst = []
for x in resultantList:
    protein_sequences = get_protein_fasta_sequence(x)
    lst.append(protein_sequences)


testLst = []
for x in test_0_values:
    protein_sequences = get_protein_fasta_sequence(x)
    testLst.append(protein_sequences)
# print(lst)

seqMatrix = []
for x in lst:
    for i in x:
        seqMatrix.append(i)
        break

testMatrix = []
for x in testLst:
    for i in x:
        testMatrix.append(i)
        break

# print(seqMatrix)
# print(len(seqMatrix))

"""Main coding

a)Mapping
"""

dicfact1 = {'A': -0.591, 'C': -1.343, 'D': 1.050, 'E': 1.375, 'F': -1.006, 'G': -0.384, 'H': 0.036, 'I': -1.29, 'K': 1.831, 'L': -
            1.019, 'M': -0.663, 'N': 0.945, 'P': 0.189, 'Q': 0.931, 'R': 1.538, 'S': -0.228, 'T': -0.032, 'V': -1.337, 'W': -0.595, 'Y': 0.260}


dicfact2 = {'A': -1.302, 'C': 0.465, 'D': 0.302, 'E': -1.453, 'F': -0.590, 'G': 1.652, 'H': -0.417, 'I': -0.547, 'K': -0.561, 'L': -
            0.987, 'M': -1.524, 'N': 0.828, 'P': 2.081, 'Q': -0.179, 'R': -0.055, 'S': 1.399, 'T': 0.326, 'V': -0.279, 'W': 0.009, 'Y': 0.830}


dicfact3 = {'A': -0.733, 'C': -0.862, 'D': -3.656, 'E': 1.477, 'F': 1.891, 'G': 1.330, 'H': -1.673, 'I': 2.131, 'K': 0.533, 'L': -
            1.505, 'M': 2.219, 'N': 1.299, 'P': -1.628, 'Q': -3.005, 'R': 1.502, 'S': -4.760, 'T': 2.213, 'V': -0.544, 'W': 0.672, 'Y': 3.097}


dicfact4 = {'A': 1.570, 'C': -1.020, 'D': -0.259, 'E': 0.113, 'F': -0.397, 'G': 1.045, 'H': -1.474, 'I': 0.393, 'K': -0.277, 'L': 1.266,
            'M': -1.005, 'N': -0.169, 'P': 0.421, 'Q': -0.503, 'R': 0.440, 'S': 0.670, 'T': 0.908, 'V': 1.242, 'W': -2.128, 'Y': -0.838}


dicfact5 = {'A': -0.146, 'C': -0.255, 'D': -3.242, 'E': -0.837, 'F': 0.412, 'G': 2.064, 'H': -0.078, 'I': 0.816, 'K': 1.648, 'L': -
            0.912, 'M': 1.212, 'N': 0.933, 'P': -1.392, 'Q': -1.853, 'R': 2.897, 'S': -2.647, 'T': 1.313, 'V': -1.262, 'W': -0.184, 'Y': 1.512}

"""b)Storing Mapping of sequences"""

matrix = list()
for a in range(len(seqMatrix)):
    s = seqMatrix[a]
    vec1 = list()
    vec2 = list()
    vec3 = list()
    vec4 = list()
    vec5 = list()
    for i in range(len(s)):
        vec1.append((dicfact1[s[i]]))
        vec2.append((dicfact2[s[i]]))
        vec3.append((dicfact3[s[i]]))
        vec4.append((dicfact4[s[i]]))
        vec5.append((dicfact5[s[i]]))

    res = []
    res.append(vec1)
    res.append(vec2)
    res.append(vec3)
    res.append(vec4)
    res.append(vec5)
    matrix.append(res)


testResultMatrix = list()
for a in range(len(testMatrix)):
    s = testMatrix[a]
    vec1 = list()
    vec2 = list()
    vec3 = list()
    vec4 = list()
    vec5 = list()
    for i in range(len(s)):
        vec1.append((dicfact1[s[i]]))
        vec2.append((dicfact2[s[i]]))
        vec3.append((dicfact3[s[i]]))
        vec4.append((dicfact4[s[i]]))
        vec5.append((dicfact5[s[i]]))

    res = []
    res.append(vec1)
    res.append(vec2)
    res.append(vec3)
    res.append(vec4)
    res.append(vec5)
    testResultMatrix.append(res)

# print(matrix)

"""C)FFT calculation
---


"""

for i in range(len(matrix)):
    for j in range(5):
        matrix[i][j] = np.fft.fft(matrix[i][j])

for i in range(len(testResultMatrix)):
    for j in range(5):
        testResultMatrix[i][j] = np.fft.fft(testResultMatrix[i][j])

"""D)Normalization of FFT sequence"""

for i in range(len(matrix)):
    for j in range(5):
        lenx = len(matrix[i][j])
        for k in range(1, lenx):
            a = np.real(matrix[i][j][k])
            b = np.imag(matrix[i][j][k])
            result = np.fabs((a**2) + (b**2))
            result = result / (math.sqrt(lenx))
            result = result/2
            matrix[i][j][k-1] = result

for i in range(len(testResultMatrix)):
    for j in range(5):
        lenx = len(testResultMatrix[i][j])
        for k in range(1, lenx):
            a = np.real(testResultMatrix[i][j][k])
            b = np.imag(testResultMatrix[i][j][k])
            result = np.fabs((a**2) + (b**2))
            result = result / (math.sqrt(lenx))
            result = result/2
            testResultMatrix[i][j][k-1] = result

"""E)Zscore Calculation"""

for i in range(len(matrix)):
    for j in range(5):
        matrix[i][j] = stats.zscore(matrix[i][j])

for i in range(len(matrix)):
    for j in range(5):
        lenx = len(matrix[i][j])
        for k in range(lenx):
            matrix[i][j][k] = np.real(matrix[i][j][k])


for i in range(len(testResultMatrix)):
    for j in range(5):
        testResultMatrix[i][j] = stats.zscore(matrix[i][j])

for i in range(len(testResultMatrix)):
    for j in range(5):
        lenx = len(testResultMatrix[i][j])
        for k in range(lenx):
            testResultMatrix[i][j][k] = np.real(testResultMatrix[i][j][k])

"""F)Max_Score Calculations"""

maxScores = list()
for i in range(len(matrix)):
    maxi = 0
    for j in range(5):
        maxi = max(maxi, max(matrix[i][j]))
    maxScores.append(abs(maxi))


testMaxScores = list()
for i in range(len(testResultMatrix)):
    maxi = 0
    for j in range(5):
        maxi = max(maxi, max(testResultMatrix[i][j]))
    testMaxScores.append(abs(maxi))

print((maxScores))

"""G)Thetha Ratio Calculation"""

thetaRationsMatrix = list()
testThethaRationMatrix = list()
for i in range(len(matrix)):
    testVector = list()
    for j in range(5):
        testVector.append(matrix[i][j])
    thetaRationsMatrix.append(compute_theta_ratio(
        testVector, len(matrix[i][j]))[3])

for i in range(len(testResultMatrix)):
    testVector = list()
    for j in range(5):
        testVector.append(matrix[i][j])
    testThethaRationMatrix.append(
        compute_theta_ratio(testVector, len(matrix[i][j]))[3])
# print((thetaRationsMatrix))
# print((resultRepeatNR))

# # Install necessary library
# # !pip install requests

# # # Import the requests module
# # import requests
# # import json

# # Function to Fetch UniProt ID from PDB ID
# def fetch_uniprot_id_from_pdb(pdb_id):
#     url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
#     response = requests.get(url)

#     if response.status_code == 200:
#         data = response.json()
#         try:
#             # Extract UniProt ID from the JSON response
#             for entity in data['rcsb_entry_container_identifiers']['polymer_entity_ids']:
#                 entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity}"
#                 entity_response = requests.get(entity_url)
#                 if entity_response.status_code == 200:
#                     entity_data = entity_response.json()
#                     for reference in entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids', []):
#                         return reference
#             return "UniProt ID not found"
#         except KeyError:
#             return "Unable to extract UniProt ID"
#     else:
#         return f"Error fetching data: HTTP {response.status_code}"

# # Replace 'YOUR_PDB_ID' with the actual PDB ID
# pdb_id = '1MYO'  # Replace with a PDB ID
# uniprot_id = fetch_uniprot_id_from_pdb(pdb_id)
# print("UniProt ID:", uniprot_id)

# # # Install necessary library
# # !pip install requests

# # # Import the requests module
# # import requests

# # Function to Fetch Data from UniProt
# def fetch_uniprot_data(uniprot_id):
#     url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
#     response = requests.get(url)

#     if response.status_code == 200:
#         return response.text
#     else:
#         return f"Error fetching data: HTTP {response.status_code}"

# # Function to Extract Repeat Information
# def extract_repeat_info(uniprot_data):
#     lines = uniprot_data.split("\n")
#     repeat_info = [line for line in lines if line.startswith("FT   REPEAT")]
#     return True if repeat_info else False

# # Fetch and Display Repeat Information
# uniprot_id = 'P62775'  # Replace with a valid UniProt ID
# uniprot_data = fetch_uniprot_data(uniprot_id)
# repeats = extract_repeat_info(uniprot_data)
# print((repeats))
# # For printing the exact repeat units INFO
# # for repeat in repeats:
# #     print(repeat)

"""H)SVM Based Model"""

z_max_values = np.array(maxScores)
theta_ratios = np.array(thetaRationsMatrix)
labels = np.array(resultRepeatNR)  # 0 for repeating, 1 for non-repeating

# Ensure the data is real (if needed)
# z_max_values = z_max_values.real
# theta_ratios = theta_ratios.real

# Combining features into a single array
X = np.column_stack((z_max_values, theta_ratios))

# Set test_size to 0 so that all data goes to the training set
test_size = 0
if test_size == 0:
    X_train, y_train = X, labels
    X_test, y_test = np.array([]), np.array([])  # Empty test set
else:
    X_train, X_test, y_train, y_test = train_test_split(X, labels, test_size=test_size, random_state=0)

# Create and fit the SVM model
model = svm.SVC(kernel='linear',probability=True)
model.fit(X_train, y_train)

# Save the model
with open('model.pkl', 'wb') as model_file:
    pickle.dump(model, model_file)

# Load the model
with open('model.pkl', 'rb') as model_file:
    svm_model = pickle.load(model_file)

print("Done")



# Predict on the test set
# predictions = model.predict(X_test)
# # out_test=model.predict(np.column_stack((theta_ratios,z_max_values)))
# # print("Done")
# # print(predictions,labels)
# # Evaluate the model
# new_z_max_values = np.array(testMaxScores)  # 7.037(R)#7.0239
# new_theta_ratios = np.array(testThethaRationMatrix)  # 2.108(R)#2.3809
#
# # Predict on the new input data
# new_predictions = model.predict(
#     np.column_stack((new_z_max_values, new_theta_ratios)))
#
# print("New Data Predictions:", new_predictions)
#
# confusion = confusion_matrix(y_test, predictions)
# sensitivity = confusion[0, 0] / (confusion[0, 0] + confusion[0, 1])
# specificity = confusion[1, 1] / (confusion[1, 0] + confusion[1, 1])
# fMeasure = 2*confusion[0, 0] / \
#     (confusion[1, 0] + confusion[0, 1] + 2*confusion[0, 0])
# mcc = matthews_corrcoef(y_test, predictions)
#
# print("Sensitivity:", sensitivity)
# print("Specificity:", specificity)
# print("Matthews Correlation Coefficient:", mcc)
# print("F-Measure:", fMeasure)
#
# # Plot the decision boundary
# h = 0.01
# x_min, x_max = z_max_values.min() - 1, z_max_values.max() + 1
# y_min, y_max = theta_ratios.min() - 1, theta_ratios.max() + 1
# xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
#
# # Ensure that the number of features in np.c_ matches the number of features used for training
# Z = model.predict(np.column_stack((xx.ravel(), yy.ravel())))
# Z = Z.reshape(xx.shape)
#
# plt.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.3)
# plt.scatter(z_max_values, theta_ratios, c=labels, cmap=plt.cm.coolwarm)
# plt.xlabel('Z-Max Values')
# plt.ylabel('Theta Ratios')
# plt.title('SVM Classification: Protein Sequences')
# plt.show()

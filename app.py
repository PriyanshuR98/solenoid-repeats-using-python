from flask import Flask, request, render_template
import pickle
import numpy as np
import math
import scipy.stats as stats
from collections import Counter
from Bio.PDB import PDBList, PDBParser, PPBuilder
import requests

app = Flask(__name__)

# Load the pre-trained SVM model
svm_model = pickle.load(open('model.pkl', 'rb'))

# Dictionary mappings
dicfact1 = {'A': -0.591, 'C': -1.343, 'D': 1.050, 'E': 1.375, 'F': -1.006, 'G': -0.384, 'H': 0.036, 'I': -1.29,
            'K': 1.831, 'L': -1.019, 'M': -0.663, 'N': 0.945, 'P': 0.189, 'Q': 0.931, 'R': 1.538, 'S': -0.228,
            'T': -0.032, 'V': -1.337, 'W': -0.595, 'Y': 0.260}
dicfact2 = {'A': -1.302, 'C': 0.465, 'D': 0.302, 'E': -1.453, 'F': -0.590, 'G': 1.652, 'H': -0.417, 'I': -0.547,
            'K': -0.561, 'L': -0.987, 'M': -1.524, 'N': 0.828, 'P': 2.081, 'Q': -0.179, 'R': -0.055, 'S': 1.399,
            'T': 0.326, 'V': -0.279, 'W': 0.009, 'Y': 0.830}
dicfact3 = {'A': -0.733, 'C': -0.862, 'D': -3.656, 'E': 1.477, 'F': 1.891, 'G': 1.330, 'H': -1.673, 'I': 2.131,
            'K': 0.533, 'L': -1.505, 'M': 2.219, 'N': 1.299, 'P': -1.628, 'Q': -3.005, 'R': 1.502, 'S': -4.760,
            'T': 2.213, 'V': -0.544, 'W': 0.672, 'Y': 3.097}
dicfact4 = {'A': 1.570, 'C': -1.020, 'D': -0.259, 'E': 0.113, 'F': -0.397, 'G': 1.045, 'H': -1.474, 'I': 0.393,
            'K': -0.277, 'L': 1.266, 'M': -1.005, 'N': -0.169, 'P': 0.421, 'Q': -0.503, 'R': 0.440, 'S': 0.670,
            'T': 0.908, 'V': 1.242, 'W': -2.128, 'Y': -0.838}
dicfact5 = {'A': -0.146, 'C': -0.255, 'D': -3.242, 'E': -0.837, 'F': 0.412, 'G': 2.064, 'H': -0.078, 'I': 0.816,
            'K': 1.648, 'L': -0.912, 'M': 1.212, 'N': 0.933, 'P': -1.392, 'Q': -1.853, 'R': 2.897, 'S': -2.647,
            'T': 1.313, 'V': -1.262, 'W': -0.184, 'Y': 1.512}


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
                    for reference in entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids',
                                                                                                          []):
                        return reference
            return "notFound"
        except KeyError:
            return "Unable to extract UniProt ID"
    else:
        return f"Error fetching data: HTTP {response.status_code}"


def get_protein_fasta_sequence(pdb_id):
    pdb_list = PDBList()
    temp_dir = "."

    pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir=temp_dir, file_format="pdb")

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

    return sequences[0]


def process_sequence(seq):
    # Mapping sequences
    sequence = get_protein_fasta_sequence(seq)
    vec1 = [dicfact1[res] for res in sequence]
    vec2 = [dicfact2[res] for res in sequence]
    vec3 = [dicfact3[res] for res in sequence]
    vec4 = [dicfact4[res] for res in sequence]
    vec5 = [dicfact5[res] for res in sequence]

    matrix = [vec1, vec2, vec3, vec4, vec5]

    # FFT calculation
    for j in range(5):
        matrix[j] = np.fft.fft(matrix[j])

    # Normalization of FFT sequence
    for j in range(5):
        lenx = len(matrix[j])
        for k in range(1, lenx):
            a = np.real(matrix[j][k])
            b = np.imag(matrix[j][k])
            result = np.fabs((a ** 2) + (b ** 2))
            result = result / (math.sqrt(lenx))
            result = result / 2
            matrix[j][k - 1] = result

    # Z-score calculation
    for j in range(5):
        matrix[j] = stats.zscore(matrix[j])

    for j in range(5):
        lenx = len(matrix[j])
        for k in range(lenx):
            matrix[j][k] = np.real(matrix[j][k])

    # Max Score Calculations
    max_zscores = max(max(matrix[j]) for j in range(5))

    # Theta Ratio Calculation
    theta_ratio = compute_theta_ratio(matrix, len(matrix[0]))[3]

    return np.real(max_zscores), theta_ratio


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/predict', methods=['POST'])
def predict():
    # Collecting form data
    pdb_id = request.form['pdb_id']

    # Processing the sequence
    try:
        max_zscores, theta_ratio = process_sequence(pdb_id)
    except Exception as e:
        return render_template('predict.html')

    # Preparing data for prediction
    features = np.array([max_zscores, theta_ratio]).reshape(1, -1)

    # Making prediction
    prediction = svm_model.predict_proba(features)[0]
    output = '{0:.{1}f}'.format(prediction[1], 2)
    print(output)
    if float(output) > 0.5:
        return render_template('predict.html',protein_id = pdb_id,repeats="Repeating")
    else:
        return render_template('predict.html',protein_id = pdb_id,repeats="Non Repeating")


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)

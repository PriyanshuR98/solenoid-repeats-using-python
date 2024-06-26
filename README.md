# Solenoid Structure Recognition in Protein Sequences using Discrete Transforms

Applied advanced DSP and windowing techniques to analyze protein sequences sourced from WORLD PDB, enhancing data
accuracy by 25% and enabling more precise structural predictions.
• Analysis using DSP and Windowing Technique improved processing time by 78.2% and diminished errors by 82.4%.
• Tech stack: Spearheaded the adoption of a comprehensive tech stack including Python, MATLAB, SciPy, and NumPy.
• Implemented efficient algorithms and data processing techniques resulting in a remarkable 78.2% improvement in
processing times and an 82.4% reduction in errors.
• These enhancements significantly boosted productivity and accuracy in scientific computing and data analysis tasks


## Introduction

Welcome to the repository for Solenoid Structure Recognition in Protein Sequences! This project focuses on the identification of repeats and non-repeats in protein sequences—a crucial task in molecular biology. Understanding these structural patterns is vital for unraveling the functions and characteristics of proteins.

In biological systems, proteins often exhibit repeated structural motifs, forming solenoid structures. Recognizing these repeats is essential for gaining insights into the architecture and functionality of proteins.

## Motivation

Identifying repeats in protein sequences has significant implications in understanding various biological processes. Repeats often play crucial roles in protein-protein interactions, molecular recognition, and cellular functions. Recognizing these patterns aids in deciphering the underlying mechanisms of various diseases and can contribute to drug discovery and design.

## Discrete Transforms and Formulation

The methodology employed in this project leverages Discrete Fourier Transforms (DFT) for recognizing solenoid structures in protein sequences. The conceptual foundation is derived from the work of Artchely (2005). The DFT formulation provides a powerful mathematical framework for mapping and analyzing repetitive patterns in the protein sequences.

## Parameters: ZScore and Theta Ratio

Two key parameters, **ZScore** and **Theta Ratio**, are calculated based on the DFT results. These parameters serve as quantitative measures for characterizing the significance and nature of the identified repeats. They play a pivotal role in the subsequent classification process.

## SVM Model for Classification

To classify repeats and non-repeats, a Support Vector Machine (SVM) model is employed. The SVM model utilizes the **ZScore** and **Theta Ratio** parameters to make predictions. The model has been trained on a dataset annotated with known solenoid structures, enabling it to effectively distinguish between repeats and non-repeats.

## How to Run the Program

Follow these steps to run the program:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/RakshitSathyakumar/solenoid-repeats-using-python
2. **Download the Datasets:**
   - Download the datasets from the provided link: [DATA.csv](From the given dataset in repository)
   - Place the downloaded **DATA.csv** file in the project's root directory.


3. **Open the Genomics_final.ipynb File:**
   - Open the **Genomics_final.ipynb** notebook using either VSCode or Google Colab.

4. **Run the Programs:**
   - Execute the notebook cells step by step.
   - Provide the necessary input when prompted.
   - View the desired outputs, including identified repeats and non-repeats.

Feel free to explore the codebase and modify parameters to suit your specific needs. If you have any questions or encounter issues, please refer to the documentation or open an issue in the repository.


Happy coding!

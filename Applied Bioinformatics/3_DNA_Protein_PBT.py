import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Task 1: Advanced Sequence handling - Recap
    # From DNA to Protein

    You're interested in comparing the RNASE1 protein sequence of a third organism with those of whales and horses. However, rather than having the protein sequence, you only possess the DNA sequence of the RNASE1 gene written in 5' to 3' direction, including information about the locations of introns.

    1. Write a function that takes a DNA sequence (written in 5' to 3' direction) as input, transcribes it into RNA, and returns the transcript in 5' to 3' direction. To do this, use a dictionary and a for loop. Consider base pair complementarity, the difference in the bases between DNA and RNA, and the direction of transcription. Call the function with your gene and save the transcript in a variable.

    2. You know that there are three introns in the transcript at the following locations: intron 1 (393 to 474), intron 2 (601 to 713), and intron 3 (853 to 1011). Remove the introns from the transcript using string slicing to get the mRNA sequence. Save the mRNA in a variable.

    3. You now want to translate the mRNA sequence into the protein sequence. First, you need to determine the position of the first start codon within the mRNA sequence. Write a function that takes an mRNA sequence as input and uses a for loop to locate the index of the first start codon. Call the function with your mRNA and save the index in a variable. Validate your result by slicing the mRNA sequence with this index.

    4. Next, write a function that takes an mRNA sequence as input and translates the mRNA into the protein sequence using a for loop, the codon_dict, and your function that locates the start codon. Consider that translation stops at a stop codon represented by "*" in the codon_dict. The "*" should not be part of the protein sequence. Call the function with your mRNA sequence and save the protein sequence in a variable.

    5. Remove the first Methionine (M) from the beginning of the protein sequence using string slicing to obtain the mature protein sequence.

    6. Find out which organism the RNASE1 protein sequence belongs to. Protein BLAST (https://blast.ncbi.nlm.nih.gov) is a fast sequence search tool that identifies the protein sequences most similar to your query sequence. A hit with 100% sequence coverage and 100% sequence identity is identical to your sequence. Go to the BLAST website, choose "Protein BLAST", paste your protein sequence into the search field, and click "BLAST". Can you find out which organism the sequence belongs to? | OR use their API https://ncbi.github.io/blast-cloud/dev/api.html

    7. Compare your protein sequence with those of whales and horses by counting the perfect matches between the sequences. According to the number of perfect matches, which of the two sequences is closer to your protein?

    8. Introduce a shift of one into your protein sequence by adding a "-" to the beginning of the sequence. Count the perfect matches again using your shifted sequence. How does it change the result?
    """)
    return


@app.cell
def _():
    gene = "CGGGAGTCTCGCCACACCGGGCTAGTATTAACTGAAGGGGTAAATAAAAACGGGCCCCTTGAGTACTGAT"\
           "TAGGGCGAAGGGCCTAATGCGCGTGGCTAAACTTCGGCACGAATAATTGGGTGGACGACTGATGGAACGC"\
           "GTTCCTTCTACGCCACCTGCTAGCGGACGTGAGAGAGGCGAAAAATTCTCTTAAAGGCTATACATATGCA"\
           "TCGAAGTGCACTGGCACGTACTGGCCTTCGCAGGCCACTATAATTTGCTTGTTCAGGTTACGGCTTTCGA"\
           "TGGCGACACCTTTCTAAAGGAGAGTAAAAGTGGAGGCCGACTCGATGGACAACAGATAGGTCAGACGCCG"\
           "GGGCGTTTTGGAACCCGGAAAGTTAAATGCGAGAGAGAGGCTTGGCATAGGGTAATGGCATCAACGTAAC"\
           "TACGGAGCTGTGGTTTCATACTGACAGTTCGGGTACTTTGAAGCCCCTGTTTGACGGCAGTTTGTAATGC"\
           "TCAATCTGGAGTTAGATTTATAGCAATTTGTACGTCCATTTTTACAAGTTACGTTCTCCTGATGACACAC"\
           "GGCATCAACTGCCTCACTCGCGACAGGTATCAGCCCTGTTGTGGGACGGAGCATTCAGGGCTATAACGGT"\
           "TATTTACAAGCATAGATGACGGGTGGTGTATGAAACTGCCTGTACTATAGTATACAGATTTCGGTTCATG"\
           "GATGAAAGTATTCAGGGGTTTGCAGCGCCCGGATGTCATGTCGCGGGCCTTCATCATTAGGTTGCAATAA"\
           "TTAGAGCTAGATGCAGTGGAGTGTTCAGTGTCCATATGCCCGGTGGCAAGGGTGGTGGTCCCGAATCGAG"\
           "AATCTGAATATCATCTAGTCTCGACGGGCTAGCGGGTGAGCGCTGTGTCTCTGTCTTTGGAATTTCTCTG"\
           "CGGGGGTTTCCATACTCCGGTCCTCTACCTCCGACGTCACTCGACACTAGTGTGGTGACACGAGGTCGGA"\
           "CTCACTGTCTCGTTCTGGGATAGAGTTTTTTTTTTTTTTTTTTCTTTTCGAGGACTCCACCTCTGCGGTT"\
           "GAGAGAGATCGAGCGATCACCCAACGTCCTCCACGAATGCGTACAAACAAAGAAACGACGGCAGAAGGTC"\
           "AACGAAATAGACAAGTGAACACGGGACTGAAAGTTGAGACAGAGGAAGGAGAAGGATGTCGTGAGGGGAC"\
           "GGGAGTTGTTCTACAAAACGGTTGACCGGTTCTGGACGGGACACGTCGACACCCAACTAAGGTGTGGGGG"\
           "CGGGCCGTGGGCGCAGGCGCGG"

    codon_dict = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
                  'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                  'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                  'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                  'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                  'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                  'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                  'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                  'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
                  'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                  'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                  'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                  'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
                  'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                  'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                  'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    Horse_RNASE1 = 'KESPAMKFERQHMDSGSTSSSNPTYCNQMMKRRNMTQGWCKPVNTFVHEPLADVQAICLQKNITCKNGQSNCYQSSSSMHITDCRLTSGSKYPNCAYQTSQKERHIIVACEGNPYVPVHFDASVEVST'

    Whale_RNASE1 = 'RESPAMKFQRQHMDSGNSPGNNPNYCNQMMMRRKMTQGRCKPVNTFVHESLEDVKAVCSQKNVLCKNGRTNCYESNSTMHITDCRQTGSSKYPNCAYKTSQKEKHIIVACEGNPYVPVHFDNSV'
    return Horse_RNASE1, Whale_RNASE1, codon_dict, gene


@app.function
# 1. Transcription
def transcribe(dna_sequence):
    # lookup for the transcription
    transcription_dict = {'G':'C', 'C':'G', 'T':'A', 'A':'U'}

    # return reversed transcript to get 5' to 3' direction
    return ''.join([transcription_dict[el] for el in dna_sequence])[::-1]


@app.cell
def _(gene):
    # 2. Remove introns
    transcribed_dna = transcribe(dna_sequence=gene)
    mRNA = transcribed_dna[:393] + transcribed_dna[474:601] + transcribed_dna[713:853] + transcribed_dna[1011:]
    print(f"mRNA: {mRNA}")
    return (mRNA,)


@app.cell
def _(mRNA):
    # 3. Find start codon
    first_codon_pos = int(mRNA.find("AUG"))
    print(f"First codon position: {first_codon_pos} | First codon: {mRNA[first_codon_pos:first_codon_pos+3]}")
    return (first_codon_pos,)


@app.cell
def _(codon_dict, first_codon_pos):
    # 4. Translate to mRNA
    def translate(mRNA_sequence):
        mRNA_tranlsation = ''
        for sequence_idx in range(first_codon_pos, len(mRNA_sequence)-4, 3):
            sub_sequence = mRNA_sequence[sequence_idx:sequence_idx+3]
            codon = codon_dict[sub_sequence]
            if codon == '*':
                break
            mRNA_tranlsation+=codon
        return mRNA_tranlsation
    return (translate,)


@app.cell
def _(mRNA, translate):
    # 5. Remove the first Methionine
    protein_sequence = translate(mRNA)[1:]
    return (protein_sequence,)


@app.cell
def _(protein_sequence):
    # 6. Identify organism
    print(f"Protein sequence: {protein_sequence} | 100% match for uncharacterized protein ASPSYDRAFT_205904 [Aspergillus sydowii CBS 593.65]")
    return


@app.function
# 7. Compare to horse and whale
# use perfect match function
def compute_perfect_matches(protein_sequence, Horse_RNASE1, Whale_RNASE1):
    matches_horse = 0
    matches_whale = 0
    for el_seq, el_horse, el_whale in zip(protein_sequence, Horse_RNASE1, Whale_RNASE1):
        matches_horse += el_seq==el_horse
        matches_whale += el_seq==el_whale
    normalized_matches_horse = matches_horse / min(len(protein_sequence), len(Horse_RNASE1))
    normalized_matches_whale = matches_whale / min(len(protein_sequence), len(Whale_RNASE1))

    print(f"Perfect matches scores for {protein_sequence}\nmRNA - Horse_RNASE1: {normalized_matches_horse}\nmRNA - Whale_RNASE1: {normalized_matches_whale}")
    return normalized_matches_horse, normalized_matches_whale


@app.cell
def _(Horse_RNASE1, Whale_RNASE1, protein_sequence):
    result_1 = compute_perfect_matches(protein_sequence, Horse_RNASE1, Whale_RNASE1)
    return


@app.cell
def _(Horse_RNASE1, Whale_RNASE1, protein_sequence):
    # 8. Shift the sequence and compare again
    # added "blank" element at the beginning as shift
    result_2 = compute_perfect_matches('-' + protein_sequence, Horse_RNASE1, Whale_RNASE1)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Task 2: PDB and Python

    In the lecture, it was shown how to retrieve a FASTA file from UniProt using the requests library. In the first part of this exercise, we would like to work out a Pipeline that retrieves a PDB file for a given ID, extracts the sequence data, and finally converts the three-letter sequence into a one-letter representation. For this task, you can use the PDB ID 1LH1 as an example. Remember that using google is not a shame.


    ---

    1. Retrieve a PDB file. Therefore, you first need to import the [request library](https://pypi.org/project/requests/). The URL for PDB will be https://files.rcsb.org/download/PDB_ID.pdb in which you have to replace the "PDB_ID" with an actual ID. Figure out what an [F-String](https://realpython.com/python-f-strings/) is and how it can be used here. With requests.get(URL) you can retrieve the content of the URL. Write a function that takes a single PDB id and returns the content of the file in text format.

    2. The [PDB file format](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)) is highly standardized and therefore well suited to be computationally exploited. Each line starts with a code word telling what kind of information will follow after. We want to get sequence information found after "SEQRES". Write a function that iterates (line-wise) over the textual content of the PDB file you retrieved, checks if a line starts with SEQRES, and returns the full three-letter encoded sequence. Hint: Use .split("\n") in order to iterate over the individual lines of the file.

    3. Now that you have the sequence you want to convert it into the [one-letter code](https://www.cup.uni-muenchen.de/ch/compchem/tink/as.html). Therefore, you first need to split the sequence into triplets and afterward convert each triplet into the one-letter representation. Hint: You may want to create a lookup dictionary.

    4. Write a function that creates a [fasta file](https://en.wikipedia.org/wiki/FASTA_format) from the sequence. It should take as input the sequence, a header string, and an optional comment string. Remember that each line in a FASTA file should contain at most 80 characters (does not apply for the header). Make sure that this rule is fulfilled.

    5. Get pyMol running in colab to visualize your pdb file.
    """)
    return


@app.cell
def _():
    # Get a PDB file with requests
    import requests

    def fetch_pdb_file(pdb_id):
        """
        Fetches a PDB file by its ID using the RCSB PDB API.

        Args:
        - pdb_id (str): The ID of the PDB file to fetch.

        Returns:
        - pdb_data (str): The content of the PDB file.
        """

        # URL for fetching PDB file by ID --> Should adapt the default url by the
        # pdb_id
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

        # Send GET request to fetch the PDB file
        response = requests.get(url)

        # return the response in text format
        return response.text
    return (fetch_pdb_file,)


@app.cell
def _(fetch_pdb_file):
    pdb_id = "1LH1"  # Example PDB ID
    pdb_data = fetch_pdb_file(pdb_id)

    print(pdb_data.split('\n')[:15])
    return pdb_data, pdb_id


@app.function
# Get the sequence from a PDB file

def extract_sequence_from_pdb(pdb_data):
    """
    Extracts the sequence information from a PDB file content.

    Args:
    - pdb_data (str): The content of the PDB file.

    Returns:
    - sequence (str): The sequence information extracted from the PDB file.
    """

    # Initialize an empty string to store the sequence
    sequence = ""

    # Split the PDB data by lines
    lines = pdb_data.split("\n")

    # Iterate through each line in the PDB data
    for line in lines:
        if "SEQRES" in line:
            sequence += ''.join([el for el in line.split(' ') if el][4:])

    return sequence


@app.cell
def _(pdb_data):
    sequence = extract_sequence_from_pdb(pdb_data)
    print(sequence)
    return (sequence,)


@app.function
# Convert to FASTA

def three_to_fasta(sequence, header="Sequence", comment=None):
    # Look up one letter code
    aa_codes = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    # Splitting the sequence into triplets
    one_letter_rep = ''
    # 80 char max line chunk for fasta format
    chunk_size = 64

    # Iterate through the sequence with a step size of 3
    for i in range(0, len(sequence), 3):
        one_letter_rep += aa_codes[sequence[i:i+3]]

    # Construct the Fasta header
    # optional add comment if provided
    fasta_header = f'> {header} - {comment}\n'

    # insert line breaks to limit each line to 80 characters
    sequence_chunks = [sequence[i:i + chunk_size] for i in range(0, len(sequence), chunk_size)]

    # Construct the Fasta format string combining sequence and header
    fasta_str = fasta_header + ''.join([chunk + '\n' for chunk in sequence_chunks])

    return fasta_str


@app.cell
def _(pdb_id, sequence):
    # write to fasta file
    fasta_sequence = three_to_fasta(sequence, pdb_id)
    print(fasta_sequence)
    with open(f'data/{pdb_id}.fasta', 'w') as f:
        f.write(fasta_sequence)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Visualize the .fasta file using the PyMol library
    """)
    return


@app.cell
def _():
    import py3Dmol
    return (py3Dmol,)


@app.cell
def _(mo, pdb_data, py3Dmol):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data)
    view.zoomTo({'chain':'A','resi':83})
    view.setBackgroundColor('white')
    view.setStyle({'chain':'A'},{'cartoon': {'color':'spectrum'}})
    view.addSurface(py3Dmol.VDW,{'opacity':0.5,'color':'lightblue'}, {'chain':'A'})
    view.zoomTo()

    # Use marimo iframe to display the interactive 3D visualization
    mo.iframe(view.write_html(fullpage=True), height=620)
    return


if __name__ == "__main__":
    app.run()

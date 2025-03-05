from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(sequence):
    """
    Generate the reverse complement of a DNA sequence.
    
    :param sequence: DNA sequence string
    :return: Reverse complement sequence
    """
    complement_map = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(complement_map)[::-1]


def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a list of sequences with their IDs.
    
    :param fasta_file: Path to FASTA file
    :return: List of tuples (sequence_id, sequence)
    """
    sequences = []
    
    try:
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append((record.id, str(record.seq)))
    except Exception as e:
        raise ValueError(f"Error parsing FASTA file: {e}")
    
    if not sequences:
        raise ValueError("No sequences found in the FASTA file")
    
    return sequences


def calculate_gc_content(sequence):
    """
    Calculate the GC content of a DNA sequence.
    
    :param sequence: DNA sequence string
    :return: GC content as a percentage
    """
    sequence = sequence.upper()
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_count = g_count + c_count
    
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0


def find_repeats(sequence, min_repeat_length=6):
    """
    Find repeating sequences in DNA.
    
    :param sequence: DNA sequence string
    :param min_repeat_length: Minimum length of repeats to find
    :return: List of tuples (repeat_sequence, positions)
    """
    repeats = []
    sequence = sequence.upper()
    
    for i in range(len(sequence) - min_repeat_length + 1):
        subseq = sequence[i:i + min_repeat_length]
        # Skip if this subsequence contains non-ACGT characters
        if any(base not in "ACGT" for base in subseq):
            continue
            
        # Find all occurrences
        positions = []
        pos = sequence.find(subseq, 0)
        while pos != -1:
            positions.append(pos)
            pos = sequence.find(subseq, pos + 1)
            
        if len(positions) > 1:
            repeats.append((subseq, positions))
    
    # Remove redundant entries (already counted in longer repeats)
    unique_repeats = []
    for repeat, positions in repeats:
        is_unique = True
        for unique_repeat, unique_positions in unique_repeats:
            # Check if this repeat is part of an already-found longer repeat
            if repeat in unique_repeat and all(p in unique_positions for p in positions):
                is_unique = False
                break
        if is_unique:
            unique_repeats.append((repeat, positions))
    
    return unique_repeats


def extract_region(sequence, start, end):
    """
    Extract a region from a DNA sequence.
    
    :param sequence: DNA sequence string
    :param start: Start position (0-based, inclusive)
    :param end: End position (0-based, exclusive)
    :return: Extracted sub-sequence
    """
    if start < 0 or end > len(sequence) or start >= end:
        raise ValueError("Invalid region coordinates")
    
    return sequence[start:end]


def translate_dna(sequence, table=1):
    """
    Translate a DNA sequence to protein using the specified genetic code.
    
    :param sequence: DNA sequence string
    :param table: Translation table number (default: 1, standard genetic code)
    :return: Protein sequence
    """
    # Ensure the sequence length is divisible by 3
    remainder = len(sequence) % 3
    if remainder != 0:
        sequence = sequence[:-remainder]
    
    # Translate using Biopython
    protein = str(Seq(sequence).translate(table=table))
    
    return protein

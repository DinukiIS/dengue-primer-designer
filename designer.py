import os
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt

def load_and_find_conserved_regions(fasta_file, window_size=20):
    """Parses Dengue genomes and finds 20bp sequences present across serotypes."""
    if not os.path.exists(fasta_file):
        print(f"Error: File '{fasta_file}' not found in this folder.")
        return []

    records = list(SeqIO.parse(fasta_file, "fasta"))
    total_sequences = len(records)
    print(f"Loaded {total_sequences} Dengue sequences from your merged file.")
    
    if total_sequences == 0:
        print("No sequences found in the file.")
        return []

    # Use the first sequence as our base template reference layout
    base_seq = str(records[0].seq).upper()
    potential_primers = []

    print("Scanning genome via sliding-window conservation approach...")
    
    # Sliding window: scan across the entire genome 20 bases at a time
    for i in range(len(base_seq) - window_size + 1):
        window = base_seq[i:i+window_size]

        # Filter out windows with ambiguous sequencing quality 'N' bases
        if "N" in window:
            continue

        # Count how many total sequences contain this exact 20bp window
        match_count = 0
        for record in records:
            if window in str(record.seq).upper():
                match_count += 1

        # Calculate conservation percentage across the dataset
        conservation_rate = (match_count / total_sequences) * 100

        # If it is conserved in at least 50% of our multi-serotype pool, save it
        if conservation_rate >= 50:
            potential_primers.append((window, i, conservation_rate))

    print(f"Found {len(potential_primers)} conserved {window_size}bp regional targets.")
    return potential_primers

def filter_primers(candidate_list):
    """Filters candidate sequences based on strict molecular biology rules."""
    print("Applying biochemical filters (GC Content 40-60% & Tm 52-60°C)...")
    viable_primers = []

    for seq, pos, rate in candidate_list:
        # 1. Calculate GC Content (Ideal for diagnostics is 40% - 60%)
        gc_count = seq.count("G") + seq.count("C")
        gc_content = (gc_count / len(seq)) * 100

        if 40 <= gc_content <= 60:
            # 2. Calculate Melting Temperature (Tm) using Wallace formula
            primer_tm = mt.Tm_Wallace(seq)

            # Ideal Tm for standard diagnostic RT-PCR is between 52C and 60C
            if 52 <= primer_tm <= 60:
                viable_primers.append({
                    "sequence": seq,
                    "position": pos,
                    "conservation": f"{rate:.0f}%",
                    "GC_content": f"{gc_content:.1f}%",
                    "Tm": f"{primer_tm:.1f}°C"
                })
    return viable_primers

if __name__ == "__main__":
    fasta_input = "dengue_sequences.fasta"
    
    # Run the design pipeline
    candidates = load_and_find_conserved_regions(fasta_input, window_size=20)
    final_primers = filter_primers(candidates)

    print("\n" + "="*65)
    print("          TOP DESIGNED DIAGNOSTIC PRIMER CANDIDATES          ")
    print("="*65)
    if final_primers:
        # Sort candidates to prioritize highest conservation rates
        final_primers.sort(key=lambda x: x['conservation'], reverse=True)
        
        for idx, p in enumerate(final_primers[:5]):  # Display the top 5 results
            print(f"Candidate {idx+1} | Pos: {p['position']} | Conserved: {p['conservation']} | Seq: {p['sequence']} | GC: {p['GC_content']} | Tm: {p['Tm']}")
    else:
        print("No candidates matched your biochemical constraints.")
    print("="*65)
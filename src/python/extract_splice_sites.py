def extract_sequences(event, genome_file, flank_length=5):
    parts = event.split(';')
    #print(parts)
    chromosome,coords,strand = parts[1].split(':')[1], parts[1].split(':')[2], parts[1].split(':')[-1]
   # chromosome, strand
 
    #print(chromosome)
    #print(coords)
    #print(strand)
    five_prime, three_prime = map(int, coords.split('-'))

    # Adjust start and end coordinates based on strand orientation
    if strand == '-':
        five_prime, three_prime = three_prime, five_prime
    #print(five_prime)
    #print(three_prime)
    
    # Read the reference genome
    with open(genome_file, 'r') as f:
        genome = f.read().split('>')[1:]

    # Find the chromosome sequence
    chromosome_seq = ''
    for entry in genome:
        if entry.startswith(chromosome):
            chromosome_seq = ''.join(entry.split('\n')[1:])
            break

    if not chromosome_seq:
        print(f"Chromosome {chromosome} not found in the reference genome.")
        return None

    # Extract the sequences with flanking regions
    five_prime_seq = chromosome_seq[five_prime - 1 - flank_length:five_prime - 1 + flank_length]
    three_prime_seq = chromosome_seq[three_prime - 1 - flank_length:three_prime - 1 + flank_length]

    return five_prime_seq, three_prime_seq

# File paths
events_file = "/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/alternative_splicing/alternative_splicing_A5_strict.ioe"
genome_file = "../../reference/fasta/Potra02_genome_hardmasked.fasta"

with open(events_file, 'r') as f:
    events = f.readlines()[1:]  # Skip the header

# Create and write sequences to a FASTA file
output_file = "A5_sequences.fasta"
with open(output_file, 'w') as out_f:
    for event in events:
        event_id = event.split('\t')[2]
        five_prime_seq, three_prime_seq = extract_sequences(event_id, genome_file)
        print("5' "+five_prime_seq)
        print("3' "+three_prime_seq)
        if five_prime_seq and three_prime_seq:
            # Write five prime sequence to FASTA
            out_f.write(f">{event_id}:five_prime\n{five_prime_seq}\n")
            # Write three prime sequence to FASTA
            out_f.write(f">{event_id}:three_prime\n{three_prime_seq}\n")
out_f.close()
print(f"Sequences saved to {output_file}.")


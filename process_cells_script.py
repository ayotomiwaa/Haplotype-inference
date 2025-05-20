import os
import subprocess
import glob
import sys

def run_simple_alignment(input_dir, reference_genome, output_dir, threads=48):
    """
    Simple, direct alignment of reads for each cell with clear error reporting
    Using samtools for barcode addition instead of pysam
    
    Parameters:
    - input_dir: Directory containing CNAsim output FASTQ files
    - reference_genome: Path to reference genome FASTA
    - output_dir: Directory to save the output
    - threads: Number of threads to use
    """
    # Create output directory
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory created: {output_dir}")
    except Exception as e:
        print(f"ERROR: Failed to create output directory: {e}")
        return None

    # Find all FASTQ files
    try:
        fastq_files = glob.glob(os.path.join(input_dir, "cell*.read1.fastq.gz"))
        if not fastq_files:
            print(f"ERROR: No FASTQ files found in {input_dir}")
            return None
        print(f"Found {len(fastq_files)} cell files")
    except Exception as e:
        print(f"ERROR: Failed to list FASTQ files: {e}")
        return None

    # Process each cell
    merged_bam = os.path.join(output_dir, "combined_cells.bam")
    individual_bams = []
    
    # First check if BWA index exists
    try:
        if not os.path.exists(reference_genome + ".bwt"):
            print(f"Creating BWA index for {reference_genome}...")
            subprocess.run(["bwa", "index", reference_genome], check=True)
    except Exception as e:
        print(f"ERROR: Failed to create BWA index: {e}")
        return None
    
    # Process each cell
    for fastq_file in fastq_files:
        try:
            # Get cell ID
            cell_id = os.path.basename(fastq_file).split('.')[0]
            read1_file = fastq_file
            read2_file = fastq_file.replace("read1", "read2")
            
            if not os.path.exists(read2_file):
                print(f"ERROR: Read2 file not found for {cell_id}: {read2_file}")
                continue
            
            print(f"Processing cell: {cell_id}")
            
            # Output files
            cell_dir = os.path.join(output_dir, cell_id)
            os.makedirs(cell_dir, exist_ok=True)
            
            sam_file = os.path.join(cell_dir, f"{cell_id}.sam")
            temp_bam = os.path.join(cell_dir, f"{cell_id}.temp.bam")
            sorted_bam = os.path.join(cell_dir, f"{cell_id}.sorted.bam")
            barcoded_bam = os.path.join(cell_dir, f"{cell_id}.barcoded.bam")
            
            # Skip if already processed
            if os.path.exists(barcoded_bam) and os.path.getsize(barcoded_bam) > 0:
                print(f"Cell {cell_id} already processed")
                individual_bams.append(barcoded_bam)
                continue
            
            # Step 1: Align reads
            print(f"Aligning reads for {cell_id}...")
            bwa_cmd = [
                "bwa", "mem",
                "-t", str(threads),
                reference_genome,
                read1_file,
                read2_file
            ]
            
            with open(sam_file, 'w') as f:
                subprocess.run(bwa_cmd, stdout=f, check=True)
            
            # Step 2: Convert to BAM
            print(f"Converting SAM to BAM for {cell_id}...")
            subprocess.run([
                "samtools", "view", 
                "-b", 
                "-o", temp_bam,
                sam_file
            ], check=True)
            
            # Step 3: Sort BAM
            print(f"Sorting BAM for {cell_id}...")
            subprocess.run([
                "samtools", "sort",
                "-@", str(threads),
                "-o", sorted_bam,
                temp_bam
            ], check=True)
            
            # Step 4: Add barcodes using samtools directly
            print(f"Adding barcodes for {cell_id}...")
            
            # Creating a 12-character barcode from cell ID
            # Extract numeric part and pad with zeros
            numeric_part = ''.join([c for c in cell_id if c.isdigit()])
            barcode = numeric_part.zfill(12)[:12]
            
            # Using awk to add the CB tag
            awk_command = f"'{{if(substr($0,1,1)!=\"@\"){{print $0\"\\tCB:Z:{barcode}\"}}else{{print}}}}'"
            shell_cmd = f"samtools view -h {sorted_bam} | awk {awk_command} | samtools view -b > {barcoded_bam}"
            
            print(f"Running command: {shell_cmd}")
            subprocess.run(shell_cmd, shell=True, check=True)
            
            # Step 5: Index the BAM
            print(f"Indexing BAM for {cell_id}...")
            subprocess.run(["samtools", "index", barcoded_bam], check=True)
            
            # Add to list of processed BAMs
            individual_bams.append(barcoded_bam)
            
            # Clean up temporary files to save space
            os.remove(sam_file)
            os.remove(temp_bam)
            
            print(f"Cell {cell_id} processed successfully")
            
        except Exception as e:
            print(f"ERROR processing cell {cell_id}: {e}")
            continue
    
    # Check if any cells were processed
    if not individual_bams:
        print("ERROR: No cells were successfully processed")
        return None
    
    print(f"Successfully processed {len(individual_bams)} cells")
    
    # Merge BAM files
    try:
        if len(individual_bams) == 1:
            # Just copy the single BAM file
            import shutil
            shutil.copy(individual_bams[0], merged_bam)
            shutil.copy(individual_bams[0] + ".bai", merged_bam + ".bai")
            print(f"Single BAM file copied to {merged_bam}")
        else:
            # Merge multiple BAM files
            print(f"Merging {len(individual_bams)} BAM files...")
            merge_cmd = ["samtools", "merge", "-@", str(threads), merged_bam] + individual_bams
            subprocess.run(merge_cmd, check=True)
            
            # Index merged BAM
            print("Indexing merged BAM...")
            subprocess.run(["samtools", "index", merged_bam], check=True)
        
        print(f"Merged BAM created at {merged_bam}")
        return merged_bam
        
    except Exception as e:
        print(f"ERROR: Failed to merge BAM files: {e}")
        return None

def run_chisel(bam_file, output_dir, reference_genome, threads=48):
    """
    Run CHISEL to extract RDRs and BAFs with clear error reporting
    
    Parameters:
    - bam_file: Path to barcoded BAM file
    - output_dir: Directory to save CHISEL results
    - reference_genome: Path to reference genome FASTA
    - threads: Number of threads to use
    """
    if not bam_file or not os.path.exists(bam_file):
        print("ERROR: Valid BAM file not provided")
        return None, None
    
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        print(f"CHISEL output directory created: {output_dir}")
        
        # Run CHISEL
        print("Running CHISEL...")
        cmd = [
            "chisel",
            "--input", bam_file,
            "--output", output_dir,
            "--reference", reference_genome,
            "--threads", str(threads)
        ]
        
        subprocess.run(cmd, check=True)
        
        # Check for output files
        rdr_file = os.path.join(output_dir, "rdr", "rdr.tsv")
        baf_file = os.path.join(output_dir, "baf", "baf.tsv")
        
        if os.path.exists(rdr_file) and os.path.exists(baf_file):
            print(f"CHISEL completed successfully")
            return rdr_file, baf_file
        else:
            print("ERROR: CHISEL completed but output files not found")
            return None, None
            
    except Exception as e:
        print(f"ERROR: CHISEL execution failed: {e}")
        return None, None

# Example usage
input_dir = "/alina-data1/Ezekiel/NSF_Cancer/chisel/src/data/chr22_nochrevents" 
reference_genome = "/alina-data1/Ezekiel/NSF_Cancer/chisel/src/data/HG00171.chr22.hap1.fa"
output_dir = "/alina-data1/Ezekiel/NSF_Cancer/Haplotype-inference"
chisel_results_dir = os.path.join(output_dir, "chisel_results")

# Run the pipeline
print("Starting alignment process...")
bam_file = run_simple_alignment(input_dir, reference_genome, output_dir)

if bam_file:
    print("Starting CHISEL analysis...")
    rdr_file, baf_file = run_chisel(bam_file, chisel_results_dir, reference_genome)
    
    if rdr_file and baf_file:
        print("Pipeline completed successfully")
        print(f"RDR file: {rdr_file}")
        print(f"BAF file: {baf_file}")
    else:
        print("Pipeline failed: CHISEL analysis did not complete successfully")
else:
    print("Pipeline failed: Alignment process did not complete successfully")
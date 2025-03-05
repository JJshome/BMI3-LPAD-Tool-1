import os
import shutil
from utils.logging_utils import Timer, setup_logger
from utils.validation import validate_sequence, validate_primer_parameters, validate_constraints
from heuristic_primer_design import find_primers
from specificity import specificity_screening
from batch_process import process_input_file, save_full_scores
from final_output import merge_results

def setup_directories(base_dir="./data/"):
    """
    Set up the directory structure needed for the LPAD pipeline.
    
    :param base_dir: Base directory for data
    :return: Dictionary of directory paths
    """
    # Define directory structure
    dirs = {
        "input": os.path.join(base_dir, "input"),
        "output": os.path.join(base_dir, "output"),
        "resource": os.path.join(base_dir, "resource"),
        "intermediate": os.path.join(base_dir, "output", "Intermediate_file"),
        "final": os.path.join(base_dir, "output", "Final_score"),
        "logs": os.path.join(base_dir, "output", "logs"),
    }
    
    # Create directories if they don't exist
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
    
    return dirs


def clean_output_directories(dirs=None):
    """
    Clean output directories to prepare for a new run.
    
    :param dirs: Dictionary of directory paths (from setup_directories)
    """
    if dirs is None:
        dirs = setup_directories()
    
    # List of directories to clean
    to_clean = [dirs["intermediate"], dirs["final"]]
    
    for dir_path in to_clean:
        if os.path.exists(dir_path):
            # Remove all files but keep the directory
            for filename in os.listdir(dir_path):
                file_path = os.path.join(dir_path, filename)
                if os.path.isfile(file_path):
                    os.unlink(file_path)


def run_pipeline(input_file, ref_file, output_dir, clean=True, verbose=True):
    """
    Run the complete LPAD pipeline with validation and timing.
    
    :param input_file: Path to input FASTA file
    :param ref_file: Path to reference genome FASTA file
    :param output_dir: Output directory for results
    :param clean: Whether to clean output directories before running
    :param verbose: Whether to print verbose output
    :return: Path to the final results file
    """
    # Setup logger
    logger = setup_logger()
    logger.info(f"Starting LPAD pipeline with input: {input_file}")
    
    # Setup directories
    dirs = setup_directories()
    if clean:
        logger.info("Cleaning output directories")
        clean_output_directories(dirs)
    
    # Validate input files
    if not os.path.exists(input_file):
        logger.error(f"Input file not found: {input_file}")
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    if not os.path.exists(ref_file):
        logger.error(f"Reference file not found: {ref_file}")
        raise FileNotFoundError(f"Reference file not found: {ref_file}")
    
    # Run the pipeline with timing
    with Timer("Complete Pipeline", logger):
        # Step 1: Find primers
        with Timer("Find Primers", logger):
            logger.info("Finding primers...")
            find_primers(input_file)
        
        # Step 2: Specificity screening
        with Timer("Specificity Screening", logger):
            logger.info("Performing specificity screening...")
            specificity_screening(input_file, ref_file)
        
        # Step 3: Process input file
        with Timer("Process Input", logger):
            logger.info("Processing input file...")
            process_input_file()
        
        # Step 4: Save full scores
        with Timer("Save Scores", logger):
            logger.info("Saving scores...")
            save_full_scores()
        
        # Step 5: Merge results
        with Timer("Merge Results", logger):
            logger.info("Merging final results...")
            merge_results(output_dir=output_dir)
    
    final_output_file = os.path.join(output_dir, "final_sorted_results.csv")
    if os.path.exists(final_output_file):
        logger.info(f"Pipeline completed successfully. Results saved to: {final_output_file}")
        return final_output_file
    else:
        logger.error("Final output file not found. Pipeline may have failed.")
        return None


def run_batch_pipeline(input_dir, ref_file, output_base_dir, file_pattern="*.fasta"):
    """
    Run the LPAD pipeline on multiple input files in batch mode.
    
    :param input_dir: Directory containing input FASTA files
    :param ref_file: Path to reference genome FASTA file
    :param output_base_dir: Base directory for output results
    :param file_pattern: Pattern to match input files
    :return: Dictionary mapping input files to their result paths
    """
    import glob
    
    logger = setup_logger()
    logger.info(f"Starting batch pipeline on directory: {input_dir}")
    
    if not os.path.exists(input_dir):
        logger.error(f"Input directory not found: {input_dir}")
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    # Find input files matching the pattern
    input_files = glob.glob(os.path.join(input_dir, file_pattern))
    if not input_files:
        logger.warning(f"No files matching pattern '{file_pattern}' found in {input_dir}")
        return {}
    
    logger.info(f"Found {len(input_files)} input files")
    
    # Create the base output directory
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Process each input file
    results = {}
    for input_file in input_files:
        # Create an output directory for this file
        file_name = os.path.basename(input_file)
        file_base = os.path.splitext(file_name)[0]
        output_dir = os.path.join(output_base_dir, file_base)
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Processing file: {file_name}")
        try:
            result_path = run_pipeline(input_file, ref_file, output_dir)
            results[input_file] = result_path
        except Exception as e:
            logger.error(f"Error processing {file_name}: {e}")
            results[input_file] = None
    
    logger.info(f"Batch processing complete. Processed {len(input_files)} files.")
    return results
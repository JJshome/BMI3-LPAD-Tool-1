#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from heuristic_primer_design import find_primers
from specificity import specificity_screening
from batch_process import process_input_file, save_full_scores
from final_output import merge_results
from utils.pipeline import run_pipeline, run_batch_pipeline
from utils.logging_utils import setup_logger

# Define default file paths
DEFAULT_INPUT_FILE = './data/input/example.fasta'
DEFAULT_REF_FILE = './data/resource/hg38.fa'
DEFAULT_OUTPUT_DIR = './data/output/Final_score/'

def main():
    # Configure argument parser
    parser = argparse.ArgumentParser(
        description='LAMP Primer Auto Design (LPAD) Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add arguments
    parser.add_argument('-i', '--input_file', type=str, default=DEFAULT_INPUT_FILE,
                      help='Path to the input FASTA file')
    parser.add_argument('-r', '--ref_file', type=str, default=DEFAULT_REF_FILE,
                      help='Path to the reference FASTA file')
    parser.add_argument('-o', '--output_dir', type=str, default=DEFAULT_OUTPUT_DIR,
                      help='Directory where output will be saved')
    parser.add_argument('-b', '--batch', action='store_true',
                      help='Run in batch mode with a directory of input files')
    parser.add_argument('-p', '--pattern', type=str, default='*.fasta',
                      help='File pattern for batch mode')
    parser.add_argument('-c', '--clean', action='store_true',
                      help='Clean output directories before running')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Print verbose output')
    parser.add_argument('--version', action='version', version='LPAD Tool v1.0.0')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Setup logger
    logger = setup_logger()
    logger.info("Starting LPAD Tool")
    
    # Run in appropriate mode
    try:
        if args.batch:
            # Check if input_file is actually a directory when in batch mode
            if not os.path.isdir(args.input_file):
                parser.error(f"In batch mode, input_file must be a directory: {args.input_file}")
                sys.exit(1)
            
            logger.info(f"Running in batch mode on directory: {args.input_file}")
            results = run_batch_pipeline(
                args.input_file, args.ref_file, args.output_dir, args.pattern
            )
            
            # Report batch results
            success_count = sum(1 for result in results.values() if result is not None)
            logger.info(f"Batch processing complete. Success: {success_count}/{len(results)}")
        else:
            # Run in single file mode
            logger.info(f"Processing single file: {args.input_file}")
            result_path = run_pipeline(
                args.input_file, args.ref_file, args.output_dir, args.clean, args.verbose
            )
            
            if result_path:
                logger.info(f"Processing complete. Results saved to: {result_path}")
                print(f"The final results have been saved to {result_path}")
            else:
                logger.error("Processing failed.")
                print("Processing failed. Check logs for details.")
                sys.exit(1)
    
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error: {e}")
        print(f"An error occurred: {e}")
        sys.exit(1)

    logger.info("LPAD Tool execution completed")

# Simple mode for backward compatibility
def legacy_main(input_file, ref_file, output):
    """Legacy main function for backward compatibility."""
    find_primers(input_file)
    specificity_screening(input_file, ref_file)
    process_input_file()
    save_full_scores()
    merge_results(output_dir=output)
    print(f"The final results have been saved to {output}")

if __name__ == "__main__":
    main()

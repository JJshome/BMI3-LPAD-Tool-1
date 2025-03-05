import os
import logging
import time
from datetime import datetime

def setup_logger(log_dir='./data/output/logs/', log_level=logging.INFO):
    """
    Set up a logger for the LPAD application with rotating file handlers.
    
    :param log_dir: Directory to store log files
    :param log_level: Logging level (default: INFO)
    :return: Configured logger instance
    """
    # Create log directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)
    
    # Create a logger
    logger = logging.getLogger('LPAD')
    logger.setLevel(log_level)
    
    # Clear any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()
    
    # Create a unique log file based on current timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'LPAD_{timestamp}.log')
    
    # Create file handler that logs all messages
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    
    # Create console handler with the same log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    
    # Create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


class Timer:
    """
    A simple timer utility for measuring execution time of code blocks.
    
    Example:
        with Timer() as t:
            # Code to time
        print(f"Operation took {t.elapsed:.2f} seconds")
    """
    def __init__(self, name=None, logger=None):
        self.name = name
        self.logger = logger
        self.start_time = None
        self.elapsed = 0
        
    def __enter__(self):
        self.start_time = time.time()
        if self.logger and self.name:
            self.logger.info(f"Starting timer: {self.name}")
        return self
    
    def __exit__(self, *args):
        self.elapsed = time.time() - self.start_time
        if self.logger:
            if self.name:
                self.logger.info(f"Timer {self.name} completed in {self.elapsed:.2f} seconds")
            else:
                self.logger.info(f"Operation completed in {self.elapsed:.2f} seconds")
                
    def reset(self):
        """
        Reset the timer to measure a new operation
        """
        self.start_time = time.time()
        if self.logger and self.name:
            self.logger.info(f"Resetting timer: {self.name}")
            
    def checkpoint(self, checkpoint_name=None):
        """
        Record a checkpoint time without stopping the timer
        
        :param checkpoint_name: Optional name for the checkpoint
        :return: Time elapsed since start or last checkpoint
        """
        current_time = time.time()
        checkpoint_time = current_time - self.start_time
        
        if self.logger:
            if checkpoint_name:
                self.logger.info(f"Checkpoint '{checkpoint_name}': {checkpoint_time:.2f} seconds")
            else:
                self.logger.info(f"Checkpoint reached: {checkpoint_time:.2f} seconds")
                
        return checkpoint_time

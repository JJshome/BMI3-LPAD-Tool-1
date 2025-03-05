def validate_primer_parameters(primer_params):
    """
    Validate primer design parameters against acceptable ranges and constraints.
    
    :param primer_params: Dictionary of primer parameters to validate
    :return: Boolean indicating if parameters are valid, and error message if not valid
    """
    # Define acceptable parameter ranges
    valid_ranges = {
        "PRIMER_INTERNAL_OPT_SIZE": (18, 30),
        "PRIMER_INTERNAL_MIN_SIZE": (15, 25),
        "PRIMER_INTERNAL_MAX_SIZE": (20, 35),
        "PRIMER_INTERNAL_OPT_TM": (55.0, 70.0),
        "PRIMER_INTERNAL_MIN_TM": (45.0, 65.0),
        "PRIMER_INTERNAL_MAX_TM": (60.0, 75.0),
        "PRIMER_INTERNAL_MAX_POLY_X": (3, 6),
        "PRIMER_INTERNAL_DNA_CONC": (25, 800),
        "PRIMER_INTERNAL_SALT_MONOVALENT": (1, 200),
        "PRIMER_INTERNAL_SALT_DIVALENT": (0, 15),
        "PRIMER_INTERNAL_DNTP_CONC": (0, 5),
        "PRIMER_NUM_RETURN": (1, 10000)
    }
    
    # Check each parameter against acceptable ranges
    for param, value in primer_params.items():
        if param in valid_ranges:
            min_val, max_val = valid_ranges[param]
            if not min_val <= value <= max_val:
                return False, f"Parameter {param} with value {value} is outside acceptable range {min_val}-{max_val}"
    
    # Check logical consistency of parameters
    if primer_params.get("PRIMER_INTERNAL_MIN_SIZE", 0) > primer_params.get("PRIMER_INTERNAL_OPT_SIZE", 0):
        return False, "PRIMER_INTERNAL_MIN_SIZE cannot be greater than PRIMER_INTERNAL_OPT_SIZE"
    
    if primer_params.get("PRIMER_INTERNAL_MAX_SIZE", 100) < primer_params.get("PRIMER_INTERNAL_OPT_SIZE", 0):
        return False, "PRIMER_INTERNAL_MAX_SIZE cannot be less than PRIMER_INTERNAL_OPT_SIZE"
    
    if primer_params.get("PRIMER_INTERNAL_MIN_TM", 0) > primer_params.get("PRIMER_INTERNAL_OPT_TM", 0):
        return False, "PRIMER_INTERNAL_MIN_TM cannot be greater than PRIMER_INTERNAL_OPT_TM"
    
    if primer_params.get("PRIMER_INTERNAL_MAX_TM", 100) < primer_params.get("PRIMER_INTERNAL_OPT_TM", 0):
        return False, "PRIMER_INTERNAL_MAX_TM cannot be less than PRIMER_INTERNAL_OPT_TM"
    
    return True, "All parameters are valid"


def validate_sequence(sequence):
    """
    Validate the input DNA sequence
    
    :param sequence: DNA sequence string to validate
    :return: Boolean indicating if sequence is valid, and error message if not valid
    """
    # Check for empty sequence
    if not sequence:
        return False, "Sequence cannot be empty"
    
    # Check length
    if len(sequence) < 50:
        return False, "Sequence is too short for LAMP primer design (minimum 50 bases)"
    
    # Check for valid DNA characters
    valid_bases = set('ACGTN')
    sequence_upper = sequence.upper()
    invalid_chars = [char for char in sequence_upper if char not in valid_bases]
    
    if invalid_chars:
        unique_invalid = set(invalid_chars)
        return False, f"Sequence contains invalid characters: {', '.join(unique_invalid)}"
    
    # Check for excessive ambiguous bases (N)
    n_count = sequence_upper.count('N')
    n_percent = (n_count / len(sequence)) * 100
    
    if n_percent > 10:
        return False, f"Sequence contains excessive ambiguous bases (N): {n_percent:.1f}%"
    
    return True, "Sequence is valid"


def validate_constraints(constraints):
    """
    Validate constraints for primer selection
    
    :param constraints: Dictionary of constraint parameters
    :return: Boolean indicating if constraints are valid, and error message if not valid
    """
    # Define acceptable constraint ranges
    valid_ranges = {
        "signature_max_length": (100, 1000),
        "min_primer_spacing": (1, 50),
        "loop_min_gap": (1, 50),
        "min_inner_pair_spacing": (1, 100),
        "opt_inner_pair_spacing": (20, 80)
    }
    
    # Check each constraint against acceptable ranges
    for param, value in constraints.items():
        if param in valid_ranges:
            min_val, max_val = valid_ranges[param]
            if not min_val <= value <= max_val:
                return False, f"Constraint {param} with value {value} is outside acceptable range {min_val}-{max_val}"
    
    # Check logical consistency of constraints
    if constraints.get("min_inner_pair_spacing", 0) > constraints.get("opt_inner_pair_spacing", 0):
        return False, "min_inner_pair_spacing cannot be greater than opt_inner_pair_spacing"
    
    return True, "All constraints are valid"

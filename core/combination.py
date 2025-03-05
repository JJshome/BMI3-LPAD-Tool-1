import csv
from core.overlap import reduce_result_by_overlap

def find_best_forward_combinations(
        inner_forward, loop_forward, middle_forward, outer_forward,
        signature_max_length, min_primer_spacing, loop_min_gap,
        include_loop_primers, ideal_gap, distance_penalties,
        penalty_weights,  # [inner, loop, middle, outer]
        to_penalty_weights  # [inner_to_loop, loop_to_middle, middle_to_outer]
):
    """
    Find the best forward combinations based on penalties and spacing constraints.

    :param inner_forward: List of inner primers.
    :param loop_forward: List of loop primers.
    :param middle_forward: List of middle primers.
    :param outer_forward: List of outer primers.
    :param signature_max_length: Maximum allowed length of a signature.
    :param min_primer_spacing: Minimum spacing between primers.
    :param loop_min_gap: Minimum gap for loop primers.
    :param include_loop_primers: Whether to include loop primers.
    :param distance_penalties: List of penalties for inter-primer distances.
    :param penalty_weights: List of weights for [inner, loop, middle, outer] primers.
    :param to_penalty_weights: List of weights for [inner_to_loop, loop_to_middle, middle_to_outer] distances.
    :return: Tuple of best forward infos, penalties, and the count of forward sets.
    """
    # Unpack weights for individual penalties
    inner_penalty_weight, loop_penalty_weight, middle_penalty_weight, outer_penalty_weight = penalty_weights

    # Unpack weights for inter-primer penalties
    inner_to_loop_penalty_weight, loop_to_middle_penalty_weight, middle_to_outer_penalty_weight = to_penalty_weights

    best_forward_infos = []
    best_forward_penalties = []
    forward_set_count = 0

    for inner_index, inner_info in enumerate(inner_forward):
        inner_location = inner_info['position']
        inner_length = inner_info['length']
        inner_penalty = inner_info['penalty']

        best_set_penalty = float('inf')  # Start with a large value

        # Calculate search range for loop primers
        search_start_at = max(0, inner_location - signature_max_length + inner_length + 20)
        loop_start_at = search_start_at
        loop_end_at = max(0, inner_location - 1 - min_primer_spacing)

        # Handle loop primer placeholder if loop primers are not included
        if not include_loop_primers:
            placeholder_primer = {'position': loop_end_at + 1, 'length': 1, 'penalty': 0}
            loop_forward = [placeholder_primer]

        for loop_index, loop_info in enumerate(loop_forward):
            loop_location = loop_info['position']
            loop_length = loop_info['length']
            loop_penalty = loop_info['penalty']

            # Check loop primer position range
            if loop_location < loop_start_at and loop_length != 1:
                continue
            if loop_location > loop_end_at and loop_length != 1:
                break

            # Calculate range for middle primers
            middle_start_at = search_start_at
            middle_end_at = min(
                loop_location - loop_length - min_primer_spacing,
                inner_location - loop_min_gap - 1
            )
            middle_end_at = max(0, middle_end_at)

            inner_to_loop_distance = inner_location - (loop_location + 1)

            for middle_index, middle_info in enumerate(middle_forward):
                middle_location = middle_info['position']
                middle_length = middle_info['length']
                middle_penalty = middle_info['penalty']

                # Check middle primer position range
                if middle_location < middle_start_at:
                    continue
                if middle_location > middle_end_at:
                    break

                # Check spacing constraints for middle primers
                if (middle_location + middle_length + min_primer_spacing > loop_location - loop_length + 1) or \
                        (middle_location + middle_length + loop_min_gap > inner_location):
                    continue

                # Calculate range for outer primers
                outer_start_at = search_start_at
                outer_end_at = middle_location - 1 - min_primer_spacing

                loop_to_middle_distance = (loop_location - loop_length + 1) - (middle_location + middle_length)

                for outer_index, outer_info in enumerate(outer_forward):
                    outer_location = outer_info['position']
                    outer_length = outer_info['length']
                    outer_penalty = outer_info['penalty']

                    # Check outer primer position range
                    if outer_location < outer_start_at:
                        continue
                    if outer_location > outer_end_at:
                        break

                    # Check spacing constraints for outer primers
                    if outer_location + outer_length + min_primer_spacing > middle_location:
                        continue

                    middle_to_outer_distance = middle_location - (outer_location + outer_length)
                    inner_to_middle_distance = inner_location - (middle_location + middle_length)

                    # Calculate penalties
                    if include_loop_primers:
                        spacing_penalty = (
                                distance_penalties[abs(inner_to_loop_distance-ideal_gap)] * inner_to_loop_penalty_weight +
                                distance_penalties[abs(loop_to_middle_distance-ideal_gap)] * loop_to_middle_penalty_weight +
                                distance_penalties[abs(middle_to_outer_distance-ideal_gap)] * middle_to_outer_penalty_weight
                        )
                        primer3_penalty = (
                                inner_penalty * inner_penalty_weight +
                                loop_penalty * loop_penalty_weight +
                                middle_penalty * middle_penalty_weight +
                                outer_penalty * outer_penalty_weight
                        )
                    else:
                        spacing_penalty = (
                                distance_penalties[inner_to_middle_distance - 30] * inner_to_loop_penalty_weight +
                                distance_penalties[abs(middle_to_outer_distance-ideal_gap)] * middle_to_outer_penalty_weight
                        )
                        primer3_penalty = (
                                inner_penalty * inner_penalty_weight +
                                middle_penalty * middle_penalty_weight +
                                outer_penalty * outer_penalty_weight
                        )

                    forward_set_penalty = spacing_penalty + primer3_penalty

                    # Update best combination if penalty is lower
                    if forward_set_penalty < best_set_penalty:
                        best_set_penalty = forward_set_penalty
                        while len(best_forward_infos) <= inner_index:
                            best_forward_infos.append(None)
                        while len(best_forward_penalties) <= inner_index:
                            best_forward_penalties.append(None)
                        best_forward_infos[inner_index] = [loop_info, middle_info, outer_info]
                        best_forward_penalties[inner_index] = [spacing_penalty, primer3_penalty]

    forward_set_count = sum(1 for item in best_forward_infos if item is not None)

    return best_forward_infos, best_forward_penalties, forward_set_count
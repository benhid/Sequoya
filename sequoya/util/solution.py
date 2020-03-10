from sequoya.problem import MSA


def restore_objs(front, problem: MSA):
    for solution in front:
        for i in range(problem.number_of_objectives):
            if not problem.score_list[i].is_minimization():
                solution.objectives[i] = -1.0 * solution.objectives[i]

    return front


def get_representative_set(front):
    """ Returns three solutions from any given front: one from the middle (by sorting the front in regard to the first
    objective) and one from each extreme region.
    """
    # find extreme regions
    upper_extreme, lower_extreme = front[0], front[-1]

    # find middle
    def _obj(s):
        return s.objectives[0]
    front.sort(key=_obj)
    middle = front[len(front) // 2]

    return upper_extreme, middle, lower_extreme

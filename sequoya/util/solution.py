from sequoya.problem import MSA


def restore_objs(front, problem: MSA):
    for solution in front:
        for i in range(problem.number_of_objectives):
            if not problem.score_list[i].is_minimization():
                solution.objectives[i] = -1.0 * solution.objectives[i]

    return front

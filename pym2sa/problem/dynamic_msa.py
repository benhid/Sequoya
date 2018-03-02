from pym2sa.core.problem import MSAProblem
from pym2sa.core.solution import MSASolution


class MSA(MSAProblem):
    def __init__(self) -> None:
        super().__init__()
        self.number_of_objectives = None
        self.number_of_variables = None
        self.number_of_constraints = None
        self.score_list = []

    def evaluate(self, solution: MSASolution) -> None:
        # self.pool.map(lambda solution: Evaluator[S].evaluate_solution(solution, problem), solution_list)

        for i in range(0, solution.number_of_objectives):
            if self.score_list[i].is_minimization():
                solution.objectives[i] = self.score_list[i].compute(solution.decode_alignment_as_list_of_pairs())
            else:
                solution.objectives[i] = self.score_list[i].compute(solution.decode_alignment_as_list_of_pairs())

    def create_solution(self) -> None:
        pass

    def get_name(self) -> str:
        return "Multiple Sequence Alignment"

from jmetal.core.problem import Problem

from pym2sa.core.solution import MSASolution


class MSAProblem(Problem[MSASolution]):
    """ Class representing MSA problems """

    def evaluate(self, solution: MSASolution) -> None:
        for i in range(self.number_of_objectives):
            if self.objectives[i].is_a_minimization_objective():
                solution.objectives[i] = self.objectives[i].compute(solution, self)
            else:
                solution.objectives[i] = -1.0 * self.objectives[i].compute(solution, self)

    def create_solution(self) -> None:
        pass

from jmetal.core.problem import Problem

from pym2sa.core.solution import MSASolution


class MSAProblem(Problem[MSASolution]):
    """ Class representing MSA problems """

    def evaluate(self, solution: MSASolution) -> None:
        pass

    def create_solution(self) -> None:
        pass

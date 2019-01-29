from abc import ABCMeta, abstractmethod

from jmetal.core.problem import Problem

from sequoya.core.solution import MSASolution


class MSAProblem(Problem[MSASolution]):
    """ Class representing MSA problems """

    __metaclass__ = ABCMeta

    def __init__(self):
        super(MSAProblem, self).__init__()
        self.number_of_constraints = 0

    @abstractmethod
    def evaluate(self, solution: MSASolution) -> MSASolution:
        pass

    def get_name(self) -> str:
        return 'Multiple Sequence Alignment problem'

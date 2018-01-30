from jmetal.core.solution import Solution


class MSASolution(Solution[str]):
    """ Class representing MSA solutions """

    def __init__(self, number_of_variables: int, number_of_objectives: int, number_of_constraints: int) -> None:
        super(MSASolution, self).__init__(number_of_variables, number_of_objectives, number_of_constraints)
        self.header = [[] for x in range(self.number_of_variables)]

    def get_length_of_alignment(self) -> int:
        return len(self.variables[0])
import time
from typing import List, TypeVar

from distributed import as_completed, Client
from jmetal.core.algorithm import Algorithm
from jmetal.core.operator import Mutation, Crossover, Selection
from jmetal.core.problem import Problem
from jmetal.operator import RankingAndCrowdingDistanceSelection
from jmetal.util.termination_criterion import TerminationCriterion

from sequoya.core.solution import MSASolution

S = TypeVar('S')
R = TypeVar('R')


def reproduction(mating_population: List[MSASolution], problem, crossover_operator, mutation_operator) -> MSASolution:
    offspring_pool = []
    for parents in zip(*[iter(mating_population)] * 2):
        offspring_pool.append(crossover_operator.execute(parents))

    offspring_population = []
    for pair in offspring_pool:
        for solution in pair:
            mutated_solution = mutation_operator.execute(solution)
            offspring_population.append(mutated_solution)

    return problem.evaluate(offspring_population[0])


class DistributedNSGAII(Algorithm[S, R]):

    def __init__(self,
                 problem: Problem[MSASolution],
                 population_size: int,
                 mutation: Mutation[MSASolution],
                 crossover: Crossover[MSASolution, MSASolution],
                 selection: Selection[List[MSASolution], MSASolution],
                 termination_criterion: TerminationCriterion,
                 number_of_cores: int,
                 client: Client):
        super(DistributedNSGAII, self).__init__()
        self.problem = problem
        self.population_size = population_size
        self.mutation_operator = mutation
        self.crossover_operator = crossover
        self.selection_operator = selection

        self.termination_criterion = termination_criterion
        self.observable.register(termination_criterion)

        self.number_of_cores = number_of_cores
        self.client = client

    def create_initial_solutions(self) -> List[S]:
        return [self.problem.create_solution() for _ in range(self.number_of_cores)]

    def evaluate(self, solutions: List[S]) -> List[S]:
        return [self.client.submit(self.problem.evaluate, solution) for solution in solutions]

    def stopping_condition_is_met(self) -> bool:
        return self.termination_criterion.is_met

    def get_observable_data(self) -> dict:
        ctime = time.time() - self.start_computing_time
        return {'PROBLEM': self.problem, 'EVALUATIONS': self.evaluations, 'SOLUTIONS': self.get_result(), 'COMPUTING_TIME': ctime}

    def init_progress(self) -> None:
        self.evaluations = self.number_of_cores

        observable_data = self.get_observable_data()
        self.observable.notify_all(**observable_data)

    def step(self) -> None:
        pass

    def update_progress(self):
        observable_data = self.get_observable_data()
        self.observable.notify_all(**observable_data)

    def run(self):
        """ Execute the algorithm. """
        self.start_computing_time = time.time()

        # create initial population
        population_to_evaluate = self.create_initial_solutions()
        task_pool = as_completed(self.evaluate(population_to_evaluate))

        self.init_progress()

        auxiliar_population = []
        for future in task_pool:
            # if initial population is not full
            if len(auxiliar_population) < self.population_size:
                received_solution = future.result()
                auxiliar_population.append(received_solution)

                new_task = self.client.submit(self.problem.evaluate, self.problem.create_solution())
                task_pool.add(new_task)
            # perform an algorithm step to create a new solution to be evaluated
            else:
                offspring_population = []

                if not self.stopping_condition_is_met():
                    offspring_population.append(future.result())

                    # replacement
                    join_population = auxiliar_population + offspring_population
                    auxiliar_population = RankingAndCrowdingDistanceSelection(self.population_size).execute(join_population)

                    # selection
                    mating_population = []

                    for _ in range(2):
                        solution = self.selection_operator.execute(population_to_evaluate)
                        mating_population.append(solution)

                    # Reproduction and evaluation
                    new_task = self.client.submit(reproduction, mating_population, self.problem,
                                                  self.crossover_operator, self.mutation_operator)

                    task_pool.add(new_task)

                    # update progress
                    self.evaluations += 1
                    self.solutions = auxiliar_population

                    self.update_progress()
                else:
                    # gratefully cancel pending solution if it's not yet running
                    future.cancel()

        self.total_computing_time = time.time() - self.start_computing_time
        self.solutions = auxiliar_population

    def get_result(self) -> R:
        return self.solutions

    def get_name(self) -> str:
        return 'dNSGA-II'

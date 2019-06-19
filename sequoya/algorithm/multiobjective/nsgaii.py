import logging
import time
from typing import List, TypeVar

import dask
from distributed import as_completed, Client
from jmetal.core.algorithm import Algorithm
from jmetal.core.operator import Mutation, Crossover, Selection
from jmetal.core.problem import Problem
from jmetal.operator import RankingAndCrowdingDistanceSelection
from jmetal.util.termination_criterion import TerminationCriterion

from sequoya.core.solution import MSASolution

LOGGER = logging.getLogger('Sequoya')

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
        return self.client.map(self.problem.evaluate, solutions)

    def stopping_condition_is_met(self) -> bool:
        return self.termination_criterion.is_met

    def get_observable_data(self) -> dict:
        ctime = time.time() - self.start_computing_time

        return {'PROBLEM': self.problem,
                'EVALUATIONS': self.evaluations,
                'SOLUTIONS': self.get_result(),
                'COMPUTING_TIME': ctime}

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

        LOGGER.info(f'Creating initial set of {self.number_of_cores}-solutions')
        create_solution = dask.delayed(self.problem.create_solution)
        evaluate_solution = dask.delayed(self.problem.evaluate)

        task_pool = as_completed([], with_results=True)

        for _ in range(self.number_of_cores):
            new_solution = create_solution()
            new_evaluated_solution = evaluate_solution(new_solution)
            future = self.client.compute(new_evaluated_solution)

            task_pool.add(future)

        batches = task_pool.batches()

        LOGGER.info(f'Creating initial population of {self.population_size} individuals')

        auxiliar_population = []
        while len(auxiliar_population) < self.population_size:
            batch = next(batches)
            for _, received_solution in batch:
                auxiliar_population.append(received_solution)

                if len(auxiliar_population) < self.population_size:
                    break

            # submit as many new tasks as we collected
            for _ in batch:
                new_solution = create_solution()
                new_evaluated_solution = evaluate_solution(new_solution)
                future = self.client.compute(new_evaluated_solution)

                task_pool.add(future)

        LOGGER.info(f'Running main loop at {time.time() - self.start_computing_time}')
        self.init_progress()

        # perform an algorithm step to create a new solution to be evaluated
        while not self.stopping_condition_is_met():
            batch = next(batches)

            for _, received_solution in batch:
                offspring_population = [received_solution]

                # replacement
                join_population = auxiliar_population + offspring_population
                auxiliar_population = RankingAndCrowdingDistanceSelection(self.population_size).execute(
                    join_population)

                # selection
                mating_population = []
                for _ in range(2):
                    solution = self.selection_operator.execute(auxiliar_population)
                    mating_population.append(solution)

                # Reproduction and evaluation
                new_task = self.client.submit(reproduction, mating_population, self.problem,
                                              self.crossover_operator, self.mutation_operator)
                task_pool.add(new_task)

                # update progress
                self.evaluations += 1
                self.solutions = auxiliar_population

                self.update_progress()

                if self.stopping_condition_is_met():
                    break

        self.total_computing_time = time.time() - self.start_computing_time

        # at this point, computation is done
        for future, _ in task_pool:
            future.cancel()

    def get_result(self) -> R:
        return self.solutions

    def get_name(self) -> str:
        return 'dNSGA-II'

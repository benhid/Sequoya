import random
from typing import List, TypeVar
import logging
import time

from jmetal.core.algorithm import Algorithm
from jmetal.core.operator import Mutation, Crossover, Selection
from jmetal.core.problem import Problem
from jmetal.operator.selection import RankingAndCrowdingDistanceSelection
from dask.distributed import Client, as_completed

logger = logging.getLogger('pyM2SA')

S = TypeVar('S')
R = TypeVar(List[S])


"""
.. module:: dNSGAII
   :platform: Unix, Windows
   :synopsis: Implementation of dNSGA-II.

.. moduleauthor:: Antonio Benítez-Hidalgo <antonio.b@uma.es>
"""


def reproduction(population: List[S], problem, crossover_operator, mutation_operator) -> S:
    if len(population) > 2:
        mating_population = random.sample(set(population), 2)
    else:
        mating_population = population

    offspring = crossover_operator.execute(mating_population)[0]
    offspring = mutation_operator.execute(offspring)

    return problem.evaluate(offspring)


class dNSGAII(Algorithm[S, R]):

    def __init__(self, population_size: int, problem: Problem[S], max_evaluations: int, mutation: Mutation[S],
                 crossover: Crossover[S, S], selection: Selection[List[S], S], number_of_cores: int, client: Client):
        super().__init__()
        self.problem = problem
        self.max_evaluations = max_evaluations
        self.mutation_operator = mutation
        self.crossover_operator = crossover
        self.selection_operator = selection

        self.population = []
        self.population_size = population_size
        self.number_of_cores = number_of_cores
        self.client = client

    def update_progress(self, population):
        observable_data = {'evaluations': self.evaluations,
                           'computing time': self.get_current_computing_time(),
                           'population': population,
                           'reference_front': self.problem.reference_front}

        self.observable.notify_all(**observable_data)

    def create_initial_population(self) -> List[S]:
        population = []

        for _ in range(self.number_of_cores):
            population.append(self.problem.create_solution())

        return population

    def run(self):
        start_computing_time = time.time()

        distributed_create_solution = self.client.scatter(self.problem.create_solution(), broadcast=True)
        distributed_problem = self.client.scatter(self.problem, broadcast=True)

        logger.debug('Creating initial population')

        population = self.create_initial_population()

        futures = []
        for solution in population:
            futures.append(self.client.submit(self.problem.evaluate, solution))

        self.evaluations += len(population)
        task_pool = as_completed(futures)

        logger.debug('Running main loop')

        while self.evaluations < self.max_evaluations:
            for future in task_pool:
                # The initial population is not full
                if len(population) < self.population_size:
                    logger.debug('Creating new solution ({0} out of {1})'.format(len(population), self.population_size))

                    received_solution = future.result()
                    population.append(received_solution)

                    new_task = self.client.submit(self.problem.evaluate, distributed_create_solution)

                    task_pool.add(new_task)
                # Perform an algorithm step to create a new solution to be evaluated
                else:
                    offspring_population = []
                    if self.evaluations < self.max_evaluations:
                        offspring_population.append(future.result())

                        # Replacement
                        join_population = population + offspring_population
                        population = RankingAndCrowdingDistanceSelection(self.population_size).execute(join_population)

                        self.update_progress(population)

                        # Selection
                        mating_population = []
                        for _ in range(2):
                            solution = self.selection_operator.execute(population)
                            mating_population.append(solution)

                        # Reproduction and evaluation
                        new_task = self.client.submit(
                            reproduction, mating_population, distributed_problem,
                            self.crossover_operator, self.mutation_operator
                        )

                        task_pool.add(new_task)

                self.evaluations += 1

                if self.evaluations % 10 == 0:
                    logger.info(
                        'PopSize: {0}. Evals: {1}. Time: {2}'.format(
                            len(population), self.evaluations, self.get_current_computing_time()
                        )
                    )

        self.total_computing_time = time.time() - start_computing_time

        self.population = population
        self.client.close()

    def get_result(self) -> R:
        return self.population

    def get_name(self) -> str:
        return 'Dynamic Non-dominated Sorting Genetic Algorithm II'

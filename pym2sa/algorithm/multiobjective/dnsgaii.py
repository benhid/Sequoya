import logging
from typing import List, TypeVar, Generic
import time

from jmetal.core.algorithm import Algorithm
from jmetal.core.operator import Mutation, Crossover, Selection
from jmetal.core.problem import Problem
from jmetal.operator import RankingAndCrowdingDistanceSelection
from dask.distributed import Client, as_completed

from pym2sa.core.solution import MSASolution

logger = logging.getLogger('pyM2SA')

S = TypeVar('S')
R = TypeVar(List[S])


"""
.. module:: dNSGA-II
   :platform: Unix, Windows
   :synopsis: Implementation of dNSGA-II.

.. moduleauthor:: Antonio J. Nebro <antonio@lcc.uma.es>
"""


class dNSGAII(Algorithm[S, R]):

    def __init__(self, problem: Problem[S], max_evaluations: int, mutation: Mutation[S], crossover: Crossover[S, S],
                 selection: Selection[List[S], S], number_of_cores: int, client: Client):
        super().__init__()
        self.problem = problem
        self.population = []
        self.max_evaluations = max_evaluations
        self.mutation_operator = mutation
        self.crossover_operator = crossover
        self.selection_operator = selection

        self.number_of_cores = number_of_cores
        self.client = client

    def update_progress(self):
        self.evaluations += 1

        observable_data = {'evaluations': self.evaluations,
                           'computing time': self.get_current_computing_time(),
                           'population': self.population,
                           'reference_front': self.problem.reference_front}

        self.observable.notify_all(**observable_data)

    def create_initial_population(self) -> List[S]:
        population = []

        for _ in range(self.number_of_cores):
            population.append(self.problem.create_solution())

        return population

    def run(self):
        population = self.create_initial_population()

        start_computing_time = time.time()

        futures = []
        for solution in population:
            futures.append(self.client.submit(self.problem.evaluate, solution))

        task_pool = as_completed(futures)

        # MAIN LOOP
        while self.evaluations < self.max_evaluations:
            for future in task_pool:
                self.evaluations += 1

                # The initial population is not full
                if len(population) < self.number_of_cores:
                    received_solution = future.result()
                    population.append(received_solution)

                    new_task = self.client.submit(self.problem.evaluate, self.problem.create_solution())
                    task_pool.add(new_task)
                # Perform an algorithm step to create a new solution to be evaluated
                else:
                    offspring_population = []
                    if self.evaluations < self.max_evaluations:
                        offspring_population.append(future.result())

                        # Replacement
                        join_population = population + offspring_population
                        self.check_population(join_population)
                        population = RankingAndCrowdingDistanceSelection(self.number_of_cores).execute(join_population)

                        # Selection
                        mating_population = []
                        for i in range(self.number_of_cores):
                            solution = self.selection_operator.execute(population)
                            mating_population.append(solution)

                        # Reproduction
                        offspring_population = []
                        for i in range(0, self.number_of_cores, 2):
                            parents = []
                            for j in range(2):
                                parents.append(mating_population[i + j])

                            offspring = self.crossover_operator.execute(parents)

                            for solution in offspring:
                                self.mutation_operator.execute(solution)
                                offspring_population.append(solution)

                        solution_to_evaluate = offspring_population[0]

                        # Evaluation
                        new_task = self.client.submit(self.problem.evaluate, solution_to_evaluate)
                        task_pool.add(new_task)

                if self.evaluations % 100 == 0:
                    logger.info("PopSize: " + str(len(population)) + ". Evals: " + str(self.evaluations))

        self.total_computing_time = time.time() - start_computing_time
        self.population = population

    def check_population(self, join_population: []):
        for solution in join_population:
            if solution is None:
                raise Exception('Solution is none')

    def get_result(self) -> R:
        return self.population

    def get_name(self) -> str:
        return 'Dynamic Non-dominated Sorting Genetic Algorithm II'


def reproduction(mating_population: List[S], problem, crossover_operator, mutation_operator) -> S:
    offspring_pool = []
    for parents in zip(*[iter(mating_population)] * 2):
        offspring_pool.append(crossover_operator.execute(parents))

    offspring_population = []
    for pair in offspring_pool:
        for solution in pair:
            mutated_solution = mutation_operator.execute(solution)
            offspring_population.append(mutated_solution)

    return problem.evaluate(offspring_population[0])


class dNSGA2MSA(dNSGAII[S, R]):

    def create_initial_population(self) -> List[MSASolution]:
        return self.problem.import_instance(self.number_of_cores)

    def run(self):
        logger.info('Importing initial population')
        population = self.create_initial_population()

        self.start_computing_time = time.time()

        futures = []
        for solution in population:
            futures.append(self.client.submit(self.problem.evaluate, solution))

        self.evaluations += len(population)
        task_pool = as_completed(futures)

        logger.info('Running main loop')
        while self.evaluations < self.max_evaluations:
            for future in task_pool:
                offspring_population = []

                if self.evaluations < self.max_evaluations:
                    offspring_population.append(future.result())

                    # Replacement
                    join_population = population + offspring_population
                    self.check_population(join_population)

                    population = RankingAndCrowdingDistanceSelection(self.number_of_cores).execute(join_population)

                    # Selection
                    mating_population = []
                    for _ in range(self.number_of_cores):
                        solution = self.selection_operator.execute(population)
                        mating_population.append(solution)

                    # Reproduction and evaluation
                    new_task = self.client.submit(reproduction, mating_population, self.problem, self.crossover_operator, self.mutation_operator)

                    task_pool.add(new_task)

                logger.info("PopSize: " + str(len(population)) + ". Evals: " + str(self.evaluations))

                self.update_progress()

        self.total_computing_time = self.get_current_computing_time()
        self.population = population
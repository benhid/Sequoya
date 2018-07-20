import logging
from typing import List, TypeVar, Generic
import time

from jmetal.core.operator import Mutation, Crossover, Selection
from jmetal.core.problem import Problem
from jmetal.operator import RankingAndCrowdingDistanceSelection
from dask.distributed import Client, as_completed

from pym2sa.core.solution import MSASolution
from pym2sa.problem import MSA

logger = logging.getLogger('pyM2SA')

S = TypeVar('S')
R = TypeVar(List[S])


"""
.. module:: dNSGA-II
   :platform: Unix, Windows
   :synopsis: Implementation of dNSGA-II.

.. moduleauthor:: Antonio J. Nebro <antonio@lcc.uma.es>
"""


class dNSGAII(Generic[S, R]):

    def __init__(self,
                 problem: Problem[S],
                 population_size: int,
                 max_evaluations: int,
                 mutation: Mutation[S],
                 crossover: Crossover[S, S],
                 selection: Selection[List[S], S],
                 number_of_cores: int,
                 client: Client):
        self.problem = problem
        self.population_size = population_size
        self.max_evaluations = max_evaluations
        self.mutation_operator = mutation
        self.crossover_operator = crossover
        self.selection_operator = selection
        self.evaluations = 0
        self.population = []
        self.total_computing_time = 0

        self.number_of_cores = number_of_cores
        self.client = client

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
                if len(population) < self.population_size:
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
                        population = RankingAndCrowdingDistanceSelection(self.population_size).execute(join_population)

                        # Selection
                        mating_population = []
                        for i in range(self.population_size):
                            solution = self.selection_operator.execute(population)
                            mating_population.append(solution)

                        # Reproduction
                        offspring_population = []
                        for i in range(0, self.population_size, 2):
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


class dNSGA2MSA(dNSGAII[S, R]):

    def create_initial_population(self) -> List[MSASolution]:
        return self.problem.import_instance(self.number_of_cores)

    def run(self):
        logger.info('Importing initial population')
        population = self.create_initial_population()

        start_computing_time = time.time()

        futures = []
        for solution in population:
            futures.append(self.client.submit(self.problem.evaluate, solution))

        task_pool = as_completed(futures)

        logger.info('Running main loop')
        while self.evaluations < self.max_evaluations:
            for future in task_pool:
                self.evaluations += 1
                offspring_population = []

                if self.evaluations < self.max_evaluations:
                    offspring_population.append(future.result())

                    # Replacement
                    join_population = population + offspring_population
                    self.check_population(join_population)

                    population = RankingAndCrowdingDistanceSelection(self.population_size).execute(join_population)

                    # Selection
                    mating_population = []
                    for i in range(self.population_size):
                        solution = self.selection_operator.execute(population)
                        mating_population.append(solution)

                    # Reproduction
                    offspring_population = []
                    for i in range(0, self.population_size, 2):
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

                print(self.evaluations)

                if self.evaluations % 100 == 0:
                    logger.info("PopSize: " + str(len(population)) + ". Evals: " + str(self.evaluations))

        self.total_computing_time = time.time() - start_computing_time
        self.population = population

import time
from typing import List, Generic, TypeVar

from distributed import as_completed, Client

from jmetal.config import store
from jmetal.core.operator import Mutation, Crossover, Selection
from jmetal.core.problem import Problem
from jmetal.operator import RankingAndCrowdingDistanceSelection

S = TypeVar('S')
R = TypeVar('R')


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


class DistributedNSGAII(Generic[S, R]):

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
        self.start_computing_time = 0
        self.observable = store.default_observable

        self.population = None
        self.number_of_cores = number_of_cores
        self.client = client

    def update_progress(self, population):
        ctime = time.time() - self.start_computing_time
        observable_data = {'EVALUATIONS': self.evaluations, 'COMPUTING_TIME': ctime, 'SOLUTIONS': population,
                           'PROBLEM': self.problem}
        self.observable.notify_all(**observable_data)

    def create_initial_solutions(self) -> List[S]:
        return [self.problem.create_solution() for _ in range(self.population_size)]

    def run(self):
        population_to_evaluate = self.create_initial_solutions()

        self.start_computing_time = time.time()

        futures = []
        for solution in population_to_evaluate:
            futures.append(self.client.submit(self.problem.evaluate, solution))

        task_pool = as_completed(futures)
        population = []

        # MAIN LOOP
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
                    population = RankingAndCrowdingDistanceSelection(self.population_size).execute(join_population)

                    # Selection
                    mating_population = []
                    for _ in range(2):
                        solution = self.selection_operator.execute(population_to_evaluate)
                        mating_population.append(solution)

                    """
                    # Reproduction
                    offspring = self.crossover_operator.execute(mating_population)
                    self.mutation_operator.execute(offspring[0])

                    solution_to_evaluate = offspring[0]

                    # Evaluation
                    new_task = self.client.submit(self.problem.evaluate, solution_to_evaluate)
                    task_pool.add(new_task)
                    """

                    # Reproduction and evaluation
                    new_task = self.client.submit(reproduction, mating_population, self.problem,
                                                  self.crossover_operator, self.mutation_operator)

                    task_pool.add(new_task)
                else:
                    print("TIME: " + str(time.time() - self.start_computing_time))
                    for future in task_pool.futures:
                        future.cancel()

                    break

                self.evaluations += 1

                if self.evaluations % 10 == 0:
                    self.update_progress(population_to_evaluate)

        self.population = population

    def get_result(self) -> R:
        return self.population

    def get_name(self) -> str:
        return 'dNSGA-II'

from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import ThreadPool
from typing import TypeVar, List
from dask.distributed import Client, as_completed

from jmetal.component.evaluator import Evaluator
from jmetal.core.problem import Problem

S = TypeVar('S')


class ParallelEvaluator(Evaluator[S]):

    def __init__(self, workers=None):
        self.pool = ThreadPool(workers)

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        self.pool.map(lambda solution: Evaluator[S].evaluate_solution(solution, problem), solution_list)

        return solution_list


class MultithreadedEvaluator(Evaluator[S]):

    def __init__(self):
        self.client = Client()

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        futures = []
        for solution in solution_list:
            futures.append(self.client.submit(problem.evaluate, solution))

        evaluated_list = []
        for future in as_completed(futures):
            evaluated_list.append(future.result())

        return evaluated_list


class SubmitEvaluator(Evaluator[S]):

    def __init__(self, submit_func):
        self.submit_func = submit_func

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        futures = [self.submit_func(Evaluator.evaluate_solution, solution, problem) for solution in solution_list]
        return [f.result() for f in futures]


class ProcessPoolEvaluator(SubmitEvaluator):

    def __init__(self, processes: int=4):
        self.executor = ProcessPoolExecutor(processes)
        super(ProcessPoolEvaluator, self).__init__(self.executor.submit)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.executor.shutdown()

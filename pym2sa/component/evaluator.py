from typing import TypeVar, List
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing.pool import ThreadPool

import dask
from dask.distributed import Client, as_completed, LocalCluster

from jmetal.component.evaluator import Evaluator
from jmetal.core.problem import Problem

S = TypeVar('S')


class MapEvaluator(Evaluator[S]):

    def __init__(self, n_workers: int=4):
        self.pool = ThreadPool(n_workers)

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        self.pool.map(lambda solution: Evaluator[S].evaluate_solution(solution, problem), solution_list)

        return solution_list

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.pool.close()


class MultithreadedEvaluator(Evaluator[S]):

    def __init__(self, n_workers: int=1):
        cluster = LocalCluster(n_workers=n_workers, processes=False)
        self.client = Client(cluster)

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        futures = []
        for solution in solution_list:
            futures.append(self.client.submit(problem.evaluate, solution))

        evaluated_list = []
        for future in as_completed(futures):
            evaluated_list.append(future.result())

        return evaluated_list


class DelayedEvaluator(Evaluator[S]):

    def __init__(self):
        pass

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        delayed_results = []

        for solution in solution_list:
            delayed_results.append(dask.delayed(problem.evaluate)(solution))

        results = dask.compute(*delayed_results)

        return list(results)


class SubmitEvaluator(Evaluator[S]):

    def __init__(self, submit_func):
        self.submit_func = submit_func

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        futures = [self.submit_func(Evaluator[S].evaluate_solution, solution, problem) for solution in solution_list]
        return [f.result() for f in futures]


class ProcessPoolEvaluator(SubmitEvaluator):

    def __init__(self, processes: int=4):
        self.executor = ProcessPoolExecutor(processes)
        super(ProcessPoolEvaluator, self).__init__(self.executor.submit)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.executor.shutdown()


class ThreadPoolEvaluator(Evaluator[S]):

    def __init__(self, workers: int=4):
        self.executor = ThreadPoolExecutor(workers)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.executor.shutdown()

    def evaluate(self, solution_list: List[S], problem: Problem) -> List[S]:
        self.executor.map(lambda solution: Evaluator[S].evaluate_solution(solution, problem), solution_list)

        return solution_list

import logging

from jmetal.util.observable import Observer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpyPopulation(Observer):
    def update(self, *args, **kwargs):
        new_population = kwargs["population"]

        for i in range(0, len(new_population)):
            logger.info("{0}-individual: {1}".format(i, new_population[i].decode_alignment_as_list_of_pairs()))

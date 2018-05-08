import logging

from jmetal.util.observable import Observer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpyPopulation(Observer):
    def update(self, *args, **kwargs):
        new_population = kwargs["population"]

        for i in range(len(new_population)-1):
            logger.info("{0}-individual: objectives {1}, length {2}"
                        .format(i, new_population[i].objectives, new_population[i].get_length_of_alignment()))

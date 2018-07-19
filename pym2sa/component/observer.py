import logging
import os

logger = logging.getLogger('pyM2SA')


class WriteSequencesToFileObserver(object):

    def __init__(self, output_directory) -> None:
        self.counter = 0
        self.directory = output_directory

        if os.path.isdir(self.directory):
            logger.warning("Directory " + self.directory + " exists. Removing contents.")
            for file in os.listdir(self.directory):
                os.remove(self.directory + "/" + file)
        else:
            logger.warning("Directory " + self.directory + " does not exist. Creating it.")
            os.mkdir(self.directory)

    def update(self, *args, **kwargs):
        with open(self.directory + "/VAR." + str(self.counter), 'w') as of:
            for solution in kwargs["population"]:
                of.write(str(solution))

        self.counter += 1

import time

from pym2sa.util.graphic import ScatterPlot


def simple() -> None:
	fig = ScatterPlot('Test plot')
	fig.plot([i for i in range(0, 10)], [j for j in range(0, 10)])

	for i in range(20):
		time.sleep(5)
		fig.update()


if __name__ == '__main__':
	simple()

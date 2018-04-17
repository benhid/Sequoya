from pym2sa.util.graphic import ScatterPlot


def interactive() -> None:
    fig = ScatterPlot('Test plot')
    fig.interactive_plot()

    for i in range(0, 100):
        fig.replace_points([i])


def simple() -> None:
    fig = ScatterPlot('Test plot')
    fig.simple_plot([i for i in range(0, 10)], [j for j in range(0, 10)], name="Test plot")


if __name__ == '__main__':
    interactive()

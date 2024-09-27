import matplotlib.pyplot as plt 

def scatter_plot(
    x, 
    y, 
    x_label, 
    y_label, 
    title, 
    path_to_save_plot, 
):
    plt.scatter(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if title is not None:
        plt.title(title)
    plt.savefig(path_to_save_plot)
    plt.clf()

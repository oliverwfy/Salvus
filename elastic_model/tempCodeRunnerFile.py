


# define 2D box domain
# length in mm
x_length = 10.
y_length = 20.
x_range = np.array([0., x_length]) * 1e-3
y_range = np.array([0., y_length]) * 1e-3


domain = sn.domain.dim2.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1])
domain.plot(return_figure=True)
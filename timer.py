import timeit

def tim(n, function, parameters = None):
    s = timeit.default_timer()
    for _ in range(n):
        function(parameters)
    e = timeit.default_timer()
    return e-s

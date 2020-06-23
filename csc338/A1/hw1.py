# g(x) = x^2 +x-4. Suppose there is a small error, h in the value of x.
def g_abs_err(x, h):
    """Returns the absolute error of computing `g` at `x` if `x` is
    perturbed by a small value `h`.
    """

    return abs(-h ** 2 - 2 * x * h - h)


def g_rel_err(x, h):
    """Returns the relative error of computing `g` at `x` if `x` is
    perturbed by a small value `h`.
    """
    abs_error = g_abs_err(x, h)
    return abs_error / (x ** 2 + x - 4)


# g(x) = x^2 +x+c
def g_root_abs_err(c, h):
    """Returns the absolute error of finding the (most) positive root of `g` when
    `c` is perturbed by a small value `h`.
    """

    return abs(-np.sqrt(1 - 4 * c) + np.sqrt(1 - 4 * c - 4 * h)) / 2


def g_root_rel_err(c, h):
    """Returns the relative error of finding the (most) positive root of `g` when
    `c` is perturbed by a small value `h`.
    """
    abs_err = g_root_abs_err(c, h)
    actural = -(-1 + np.sqrt(1 - 4 * c)) / 2
    return abs_err / actural


# q2 start----------------------------------------------------------------------
# q2 start----------------------------------------------------------------------
import matplotlib.pyplot as plt
import math
import numpy as np


def f(x):
    return (x - math.sin(x)) / math.pow(x, 3)


def plot_f():
    xs = [x for x in np.arange(-3.0, 3.0, 0.05) if abs(x) > 0.05]
    ys = [f(x) for x in xs]
    plt.plot(xs, ys, 'bo')
    plt.show()


x = 0.00000001
q2_est = f(x)
first = 1 / x ** 2
second = - math.sin(x) / math.pow(x, 1)
q2_true = first + second
from decimal import *


def f2(x):
    """
    I use Maclaurin expansion (the Taylor expansion about 0)
    """

    bigger = 1 / math.pow(x, 3)
    # localcontext().prec = 100
    summ = Decimal(x * bigger)

    i = 2
    flag = 2
    while i < 1000:
        fac = Decimal(math.factorial(i))
        if flag % 2 == 0:
            summ -= Decimal(math.pow(x, i)) * Decimal(bigger) / fac
        else:
            summ += Decimal(math.pow(x, i)) * Decimal(bigger) / fac
        i += 2
        flag += 2
    return summ


# q3 start----------------------------------------------------------------------
# q3 start----------------------------------------------------------------------
xs = [0.1, 0.5, 1.0, 3.0]
q3_forward = [None, None, None, None]
q3_backward = [None, None, None, None]

true_value = [math.sin(i) for i in xs]

def function(x):
    return x - math.pow(x, 3) / math.factorial(3)
for i  in range(4):
    q3_forward[i] = true_value[i] - function(xs[i])
    q3_backward[i] = "DNE"

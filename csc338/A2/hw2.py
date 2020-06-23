# 1--------------------------------
# ---------------------------------
import math

import numpy as np


def q1():
    # beta  =3  , precision p =5 , [L,U] = [-3,3]
    BASE = 3
    PRECISION = 5
    L = -3
    U = +3
    example_float = ([1, 0, 0, 1, 0], 1, 1)
    total_num_float = (BASE ** 1 - 1) * (BASE ** PRECISION)

    print("num of  normalized-point numbers  in our system:",total_num_float)

    def is_valid_float(float_value):
        """Returns a boolean representing whether the float_value is a valid,
        normalized float in our floating point system.
        >>> is_valid_float(([1, 0, 0, 1, 0], 1, 1))
        True
        """
        (mantissa, exponent, sign) = float_value
        for i in mantissa:
            if i not in (0, 1, 2):
                return False
        if not -3 <= exponent <= 3:
            return False
        if sign != 1 and sign != -1:
            return False
        return True


    largest_negative = ([2, 2, 2, 2, 2], 3, -1)
    smallest_positive = ([1, 0, 0, 0, 0], -3, 1)
    float_32 = ([1, 0, 1, 2, 0], 3, 1)
    print("largest_negative is :",largest_negative)
    print("smallest_positive is :",smallest_positive)
    print("float_32 is :",float_32 )
    def to_float(float):
        """
        convert the number in base ten to then tuple of the form
        >>> to_float(3.111)
        ([1, 0, 0, 1, 0], 1, 1)
        """
        pass
    def to_num(float_value):
        """
        converts a tuple of the form (mantissa, exponent, sign) to a Python
        numerical value.
        >>> to_num(example_float)
        3.111111111111111
        """
        if not is_valid_float(float_value):
            return "Not a valid float number"
        (mantissa, exponent, sign) = float_value
        summ = 0
        for i in range(5):
            summ += (mantissa[i] / 3) ** i
        summ = summ * 3 ** exponent
        result = summ if sign == 1 else -summ
        return result
    def add_float(float1, float2):
        """Return a valid floating-point representation of the form (mantissa, exponent, sign)
        that is the sum of ‘float1‘ and ‘float2‘. Raises a ValueError if the result of
        the addition is not a valid float.>>> add_float(example_float, example_float)
        ([2, 0, 0, 2, 0], 1, 1])
        """
        (mantissa1, exponent1, sign1) = float1
        (mantissa2, exponent2, sign2) = float2
        # You may assume that sign1 and sign2 are positive
        assert (sign1 == 1) and (sign2 == 1)
        diff = exponent1 - exponent2
        fmantissa = []
        if diff >= 0:
            for i in range(5 - diff):
                fmantissa.append(mantissa1[i + diff] + mantissa2[i])
            for i in range(diff):
                fmantissa.insert(i, mantissa1[i])
        if diff < 0:
            for i in range(5 + diff):
                fmantissa.append(mantissa1[i] + mantissa2[i - diff])
            for i in range(abs(diff)):
                fmantissa.insert(i, mantissa2[i])
        # calculate the new mantissa
        for i in range(4, 0, -1):
            if fmantissa[i] > 2:
                fmantissa[i - 1] += 1
                fmantissa[i] = fmantissa[i] % 3
        # calculate the new exponent
        if fmantissa[0] > 2:
            final_exponent = max(exponent1, exponent2) + 1
            fmantissa[0] = fmantissa[0] % 3
        else:
            final_exponent = max(exponent1, exponent2)
        if is_valid_float((fmantissa, final_exponent, sign1)):
            return (fmantissa, final_exponent, sign1)
        else:
            print((fmantissa, final_exponent, sign1))
            raise ValueError

    f1 = ([2, 2, 2, 2, 2], 2, 1)
    f2 = ([2, 2, 2, 2, 2], 2, 1)
    # print(add_float(f1,f2))

    #----------------------------------------Q2F
    def q2f():
        a = smallest_positive  ## smalllest_postive  is too long to type in
        exmaple1f = ([0, 0, 1, 0, 0], 1, 1)  # the original number
        # print(add_float(add_float(add_float(add_float(exmaple1f, a), a), a), a))  # use (((1+a)+a)+a)+a) method
        suma = ([0, 0, 0, 0, 0], 0, 1)
        for i in range(4):
            suma = add_float(suma, a)
        print(add_float(exmaple1f, suma))  # use (1+(a+(a+(a+(+a)))) method  . they are not the same


# 2--------------------------------
# ---------------------------------
def h1(x, n):
    """Returns a list of the first n terms of the Taylor Series expansion of 1/(1-x)."""
    return [pow(x, i) for i in range(n)]


def h2(x, n):
    """Returns a list of the first n terms of the Taylor Series expansion of e^x."""
    return [pow(x, i) / math.factorial(i) for i in range(n)]


def q2():
    x = -0.01
    h11 = h1(x, 15)
    h22 = h2(x, 15)

    q2h1list = [sum(h11[0:i]) for i in range(15)]
    q2h2list = [sum(h22[0:i]) for i in range(15)]
    # print("The summation of items from 0-1 ,0-2,0-3.....  0-15 items in h1() list is :", q2h1list)
    # print("The summation of items from  0-1,0-2,0-3.....  0-15 items in h2() list is :", q2h2list)

# 3--------------------------------
# ---------------------------------
def q3():
    def z(n):
        a = pow(2.0, n) + 10.0
        b = (pow(2.0, n) + 5.0) + 5.0
        return a - b
    nonzero_zn = []
    n = 0
    while 1:
        try:
            zn = z(n)
        except OverflowError:
            break
        if zn != 0:
            nonzero_zn.append(n)
        else:
            pass
        n = n + 1
    print(nonzero_zn)  #[53, 55, 56]

        # for i in nonzero_zn:
        #     print(z(i))


if __name__ == '__main__':
    q1()
    q2()
    q3()







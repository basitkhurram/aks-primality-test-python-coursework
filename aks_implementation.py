from polynomial import Polynomial, binomial_coefficient
from math import log
from fractions import gcd

FAST_POWERING_BUFFER = 5000

def aks_prime(n):
    return aks(n) == "Prime"

def aks(n):
    if n < 2:
        return "Please don't be silly"

    if power_check(n):
        return "Composite"

    r = rfind(n)

    for a in range(2, r + 1, 1):
        hcf = gcd(a, n)
        if 1 < hcf and hcf < n:
            return "Composite"

    if n <= r:
        return "Prime"


    phi_r = totient(r)
    limit = int((phi_r**0.5) * log(n, 2))

    global fermat_generate_lhs
    if n > FAST_POWERING_BUFFER:
        fermat_generate_lhs = fermat_sm_generate_lhs
    else:
        fermat_generate_lhs = fermat_bin_generate_lhs
    for a in range(limit, 0, -1):
        if not fermat_congruence(a, n, r):
            return "Composite"

    return "Prime"


def power_check(n):
    #trying to find a (limit + 1)root
    #will result in a value less than two
    #the only positive integer less
    #than two is one
    #power_check needs to find a value
    #of two or greater

    limit = log(n, 2)
    b = 2
    while b <= limit:
        a = n**(1/b)
        if a % 1 == 0:
            return True
        b += 1
    return False

def rfind(n):
    limit = (log(n, 2))**2
    #by definition, r cannot be less than three
    r = 3
    
    while n_order_r(n, r) <= limit:
        r += 1
    return r
    
def n_order_r(n, r):
    #check for valid input
    if r < 2:
        return 0
    if n == 1:
        return 1

    k = 2
    product = pow(n, k, r)
    #to find terms that have infinite order
    old_products = set()
    #old_products = []
    #old_products = {}
    while product > 1:
        k += 1
        product = pow(n, k, r)
        if product in old_products:
            return 0
        old_products.add(product)
        #old_products.append(product)
        #old_products[product] = None
    return k

def totient(r):
    count = 0
    for n in range(1, r + 1, 1):
        if gcd(n, r) == 1:
            count += 1
    return count


def fermat_sm_generate_lhs(a, n, r):
    p = Polynomial('x + {0}'.format(a))
    p = p.__pow__(n, n, r)
    return p

def fermat_bin_generate_lhs(a, n, r):
    coefficients = [(pow(a, k, n) * binomial_coefficient(n, k)) \
                    % n for k in range(n + 1)]
    return Polynomial(coefficients)

def fermat_congruence(a, n, r):
    divisor = Polynomial('x^{0} - 1'.format(r))
    LHS = fermat_generate_lhs(a, n, r)
    RHS = Polynomial("x^{0} + {1}".format(n % r, a))
    difference = RHS - LHS
    difference = difference.polynomial_special_mod(n, r)
    return difference == Polynomial('0')

def prime(n):
    if n < 2:
        return False
    factor = 2
    limit = round(n**0.5)
    while factor < limit + 1:
        if n % factor == 0:
            return False
        factor += 1
    return True

class Polynomial:
    '''
    Polynomial class
    '''
    
    def __init__(self, coefficient_array):

        if type(coefficient_array) == str:
            #enables parsing strings as polynomials
            #dicts are used, keys are degrees, values are coefficients
            #ordered lists are used, dicts are unordered
            poly = dict()
            expression = coefficient_array
            coefficient_array = []
            exponent = ""
            term = ""
            power = False
            negative = False
            constant = True
            if expression[0] == "-":
                negative = True
                expression = expression[1:]
            powers = []
            expression = "".join(expression.split())

            for character in expression:

                if character in ["+", "-", "^", "x"]:

                    if term == "":
                        term = 1

                    if negative:
                        term = -1 * eval(term)

                    term = str(term)

                    negative = False
                
                if character == "^":
                    if term == "":
                        term = 1
                    constant = True
                    
                    continue

                if not character in ["+", "-"]:
                    
                    if power:
                        exponent += character
                        continue

                    else:
                        term += character
                        if character == "x":
                            power = True
                            constant = False
                            term = term[: -1]
                            continue

                if character in ["+", "-"]:
                    
                    if character == "-":
                        negative = True

                    if exponent == "":
                        if constant:
                            exponent = 0
                        else:
                            exponent = 1                        

                    if eval(term) % 1 == 0:
                        term = eval(term)
                    else:
                        term = float(term)
                        
                    exponent = int(exponent)

                    poly[exponent] = term
                    powers.append(exponent)                    
                    
                    power = False

                    term = ""
                    exponent = ""
                    continue            

            if term == "":
                term = 1

            if exponent == "":
                if expression[-1] == "x":
                    exponent = 1
                else:
                    exponent = 0

            if eval(term) % 1 == 0:
                term = eval(term)
            else:
                term = float(term)
            if negative:
                term = -1 * term

            exponent = int(exponent)

            poly[exponent] = term
            powers.append(exponent) 


            powers.sort(reverse = True)
            max_degree = int(powers[0])

            for degree in range(max_degree, -1, -1):
                if degree in poly:
                    coefficient_array += [poly[degree]]
                else:
                    coefficient_array += [0]   
        self._coefficients = list(coefficient_array)

        #cleans input, in case polynomial is
        #padded with redundant zeros
        while len(self._coefficients) > 0 and self._coefficients[0] == 0:
            self._coefficients.pop(0)

        self._degree = len(self._coefficients) - 1

    def get_coefficients(self):
        '''
        Returns all the coefficients of the polynomial
        '''
        return self._coefficients

    def get_coefficient(self, power):
        '''
        Returns the coefficient of the specified power
        '''
        return self.get_coefficients()[self.get_degree() - power]

    def get_degree(self):
        return self._degree

    def set_coefficient(self, power, coefficient):
        while power > self.get_degree():
            [0] + self._coefficients
        self._coefficients[self.get_degree() - power] = coefficient

    def poly_clone(self):
        return Polynomial(self.get_coefficients())

    def __repr__(self):

        '''
        Returns a string representation
        of a polynomial in the form
        ax^n + ... + bx + c,
        n is the degree of the polynomial
        a, b, and c are coefficients.
        '''

        coefficients = self.get_coefficients()
        
        output = ""

        #empty polynomial returns zero
        if self._degree == -1:
            return "0"
        
        #to keep track of the degree of the coefficient
        temp_deg = self.get_degree()

        for coefficient in coefficients:
            coefficient = eval(str(coefficient))
            #ignore zero coefficients
            if coefficient != 0:

                if coefficient < 0:
                    sign = "-"
                    coefficient = abs(coefficient)
                else:
                    sign = "+"
                    
                if coefficient == 1 and temp_deg > 0:
                    coefficient = ""
                output += " {0} {1}x^{2}".format(sign, coefficient, \
                                                 temp_deg)
                if temp_deg == 1:
                    output = output[:-2]
            temp_deg -= 1

        if output[1] == "-":
            output = "-" + output[3:]
        else:
            output = output[3:]

        if output[-3:] == "x^0":
            output = output[:-3]
        return output

    def __eq__(self, q):
        return (self - q).get_coefficients() == []

    def __add__(self, q):

        p_coefficients = self.get_coefficients()
        q_coefficients = q.get_coefficients()

        p_deg = self.get_degree()
        q_deg = q.get_degree()

        if p_deg < q_deg:
            p_coefficients, q_coefficients = q_coefficients, \
                                             p_coefficients
            p_deg, q_deg = q_deg, p_deg

        q_pseudo_deg = q_deg

        while q_pseudo_deg < p_deg:
            q_coefficients = [0] + q_coefficients
            q_pseudo_deg += 1

        p_plus_q_coefficients = list(map(lambda x,y : x + y, \
                                         p_coefficients, q_coefficients))
        p_plus_q = Polynomial(p_plus_q_coefficients)

        return p_plus_q

    def __sub__(self, q):
        neg_q = Polynomial([-1])*q
        return self.__add__(neg_q)

    def __mul__(self, q):
        if type(q) == Polynomial:
            
            p_coefficients = self.get_coefficients()
            q_coefficients = q.get_coefficients()

            p_deg = self.get_degree()
            q_deg = q.get_degree()
            pq_deg = p_deg + q_deg

            if p_deg < q_deg:
                p_coefficients, q_coefficients = q_coefficients, \
                                                 p_coefficients
                p_deg, q_deg = q_deg, p_deg

            pq_coefficients = [0 for _ in range(pq_deg + 1)]

            p_term_power = p_deg
            while p_term_power > -1:
                p_coefficient = p_coefficients[p_deg - p_term_power]
                q_term_power = q_deg
                while q_term_power > -1:
                    q_coefficient = q_coefficients[q_deg - q_term_power]
                    pq_coefficients[pq_deg - (p_term_power + q_term_power)] += \
                                           p_coefficient * q_coefficient
                    q_term_power -= 1
                p_term_power -= 1

        else:
            pq_coefficients = [coefficient*q for coefficient in \
                               self.get_coefficients()]

        return Polynomial(pq_coefficients)

    def __rmul__(self, q):
        return self.__mul__(q)

    def long_division(self, divisor):
        dividend_clone = self.poly_clone()
        divisor_clone = divisor.poly_clone()

        dividend_coefficients = dividend_clone.get_coefficients()
        divisor_coefficients = divisor_clone.get_coefficients()

        dividend = Polynomial(dividend_coefficients)
        
        if dividend.is_scalar_multiple(divisor):
            leading_dividend = dividend.get_coefficients()[0]
            leading_divisor = divisor.get_coefficients()[0]
            common_factor = leading_dividend/leading_divisor
            if common_factor % 1 == 0:
                common_factor = int(common_factor)
            quotient = Polynomial([common_factor])
            return quotient, Polynomial([0])

        deg_dividend = self.get_degree()
        deg_divisor = divisor.get_degree()
        deg_quotient = deg_dividend - deg_divisor

        quotient_coefficients = [0 for _ in range(deg_quotient + 1)]

        while deg_divisor <= deg_dividend:
            leading_dividend = dividend.get_coefficient(deg_dividend)
            leading_divisor = divisor.get_coefficient(deg_divisor)

            quotient_term = leading_dividend/leading_divisor
            if quotient_term % 1 == 0:
                quotient_term = int(quotient_term)
            term_power = deg_dividend - deg_divisor

            quotient_coefficients[deg_quotient - term_power] = \
                                               quotient_term

            quotient_term = Polynomial([quotient_term] + \
                                       [0 for _ in range(term_power)])
            subtractant = quotient_term * divisor
            dividend -= subtractant
            deg_dividend = dividend.get_degree()

        quotient = Polynomial(quotient_coefficients)
        remainder = dividend
        return quotient, remainder

    def __pow__(self, power, con_mod = None, poly_mod = None):
        if power == 0:
            return 1

        product = self
        temp = 1
        while power > 1:
            if power % 2 == 1:
                temp = (product * temp).polynomial_special_mod(con_mod, poly_mod)
            product = (product * product).polynomial_special_mod(con_mod, poly_mod)
            power //= 2
        product = (product * temp).polynomial_special_mod(con_mod, poly_mod)
        return product

    def is_scalar_multiple(self, q):
        p_deg = self.get_degree()
        q_deg = q.get_degree()
        if p_deg != q_deg:
            return False
        leading_p = float(self.get_coefficient(0))
        leading_q = float(q.get_coefficient(0))
        return leading_p * q == leading_q * self

    def polynomial_special_mod(self, n = None, r = None):
        p_degree = self.get_degree()
        p_coefficients = self.get_coefficients()
        if r != None:
            q_degree = min(r - 1, p_degree)
            q_coefficients = [coefficient for coefficient in p_coefficients[p_degree - q_degree:]]
            for k in range(0, p_degree - q_degree, 1):
                q_coefficients[-((p_degree + 1 - k) % r)] += p_coefficients[k]
        else:
            q_coefficients = self.get_coefficients()
        if type(n) == int or type(n) == float:
            q_coefficients = [coefficient % n for coefficient in q_coefficients]
        return Polynomial(q_coefficients)

computed_factorials = {}
def factorial(n):
    global computed_factorials
    if n < 2:
        return 1
    if n in computed_factorials:
        product = computed_factorials[n]
    else:
        product = 1
        for k in range(n, 1, -1):
            product *= k
    computed_factorials[n] = product
    return product

computed_binomial_coefficients = {}        
def binomial_coefficient(n, k):
    global computed_binomial_coefficients
    if (n, k) in computed_binomial_coefficients:
        return computed_binomial_coefficients[(n, k)]
    else:
        output = factorial(n) // (factorial(k) * factorial(n - k))
        computed_binomial_coefficients[(n, k)] = output
        computed_binomial_coefficients[(n, n - k)] = output
        return output

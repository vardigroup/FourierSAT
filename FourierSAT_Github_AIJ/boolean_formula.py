from collections import defaultdict
from FourierCoefficient import CardinalityFC,binom
import math
import scipy
from fractions import Fraction as F
def canonicalize(literals, k, coefs, comparator):
    if comparator == "<=":
        comparator = ">="
        coefs = [ -coefs[i] for i in range(len(coefs))]
        k = -k
    for i in range(len(coefs)):
        if coefs[i] < 0:
            coefs[i] = -coefs[i]
            literals[i] = -literals[i]
            k += coefs[i]
    return literals, k, coefs, comparator


class Formula:
    def __init__(self):
        self._iswcnf = 0
        self._variables = []
        self._clauses = []
        self._klist = []
        self._weight = []
        self._coefs = []
        self._ctype = []
        self._optsol = 0
        self._comparator = []
    def add_clause(self, literals,k, ctype, weight, coefs=None, comparator=None):
        """
        Add a new CNF clause representing the disjunction of the provided literals.

        The elements of literals should be variable id (representing the corresponding
        positive literal) or the negative of a variable id (for the negative literal).

        :param literals: An iterable of variable ids and negations of variable ids.
        :return: None
        """
        self._clauses.append(list(literals))
        self._klist.append(k)
        self._weight.append(weight)
        self._ctype.append(ctype)
        self._coefs.append(coefs)
        self._comparator.append(comparator)
    def fresh_variable(self, neg_weight, pos_weight):
        """
        Create a new weighted variable in the formula.

        :param neg_weight: Multiplicative weight on an assignment when variable is false
        :param pos_weight: Multiplicative weight on an assignment when variable is true
        :return: An id of the new variable, which can be used to construct clauses
        """
        new_var_id = len(self._variables) + 1
        self._variables[new_var_id] = (neg_weight, pos_weight)
        return new_var_id

    @property
    def clauses(self):
        return list(self._clauses)

    @property
    def variables(self):
        return self._variables

    def literal_weight(self, lit):
        """
        Returns the multiplicative weight of the provided DIMACS literal.
        """
        return self._variables[abs(lit)][1 if lit > 0 else 0]

    def write_cachet(self, filename):
        """
        Write the formula into the format expected by cachet.

        If the weights are non-probabilistic (i.e., positive and negative weights do not add up to 1),
        the weights must be normalized before they are input into cachet. This results in a normalization constant,
        which is returned.

        The proper count of the formula is then the normalization constant multiplied by the result of cachet
        when run on the output file.

        :param filename: The file to write the formula.
        :return: The normalization constant C.
        """
        with open(filename, 'w') as f:
            f.write('p cnf %d %d\n' % (len(self._variables), len(self._clauses)))

            # Write all weights, renormalizing if appropriate
            normalization_constant = 1
            for var, weight in self._variables.items():
                negative_weight = weight[0]
                positive_weight = weight[1]

                if positive_weight == negative_weight:
                    normalization_constant *= positive_weight
                    f.write('w\t%d\t%f\n' % (var, -1))  # special cachet syntax for 1-1 variables
                else:
                    # Normalize the weights, if required
                    if positive_weight + negative_weight != 1:
                        normalization_constant *= (positive_weight + negative_weight)
                        positive_weight /= (positive_weight + negative_weight)
                        f.write('w\t%d\t%f\n' % (var, positive_weight))

            # Write all clauses.
            f.writelines(['%s 0\n' % ' '.join(map(str, clause)) for clause in self._clauses])
        return normalization_constant

    def write_DIMACS(self, filename):
        """
        Write the formula into DIMACS format.

        Does not include variable weights.

        :param filename: The file to write the formula.
        :return: None
        """
        with open(filename, 'w') as f:
            f.write('p cnf %d %d\n' % (len(self._variables), len(self._clauses)))
            f.writelines(['%s 0\n' % ' '.join(map(str, clause)) for clause in self._clauses])

    @staticmethod
    def read_DIMACS(filename, include_missing_vars=False):
        """
        Read a DIMACS file into a formula.

        The file may optionally contain weights, in the cachet style i.e. lines of the form
        w [var id] [prob]
        that each indicate that the variable [var id] should have positive literal weight [prob]
        and negative literal weight 1-[prob].

        :param filename: The file to read
        :param include_missing_vars: If true, variables indicated by the DIMACS header are assigned a weight 1 1
        :return: the resulting formula
        """
        result = Formula()
        vars = defaultdict(lambda: result.fresh_variable(1, 1))

        num_vars = 0
        with open(filename, 'r') as f:
            for line in f:
                split = line.split()
                if len(split)>=2 and split[1] == '#variable=':
                    num_vars = int(line.split()[2])
                    result._variables = [i+1 for i in range(num_vars)]
                if len(line) == 0 or line[0] == 'c' or line[0] == '*':
                    continue
                if line[0] == 'p':
                    num_vars = int(line.split()[2])
                    result._variables = [i+1 for i in range(num_vars)]
                elif line[0] == 'g':
                    args = line.split()
                    k = int(args[1])
                    weight = magic_weight(num_vars,k,'g')
                    result._optsol = abs(k)
                    coefs = [1 for i in range(num_vars)]
                    if k>=0:
                        result.add_clause(range(1,num_vars+1),k,'c',weight, coefs, '>=')#num_vars)#2/(1-c))
                    else:
                        result.add_clause(range(-num_vars,0),num_vars+k,'c', weight,coefs, '>=')#num_vars)#2/(1+c))
                elif line[0] == 'x':
                    literals = map(int, [line.split()[i] for i in range(1,len(line.split())-1)])
                    literals = list(literals)
                    if len(literals) == 0:
                        continue
                    result.add_clause(literals,0,'x',magic_weight(len(literals),1,'x'))
                elif line[0] == 'n':
                    literals = map(int, [line.split()[i] for i in range(1,len(line.split())-1)])
                    literals = list(literals)
                    if len(literals) == 0:
                        continue
                    result.add_clause(literals,0,'n',magic_weight(len(literals),1,'n'))
                elif line[0] == 'd': # cardinality constraints
                    k = int(line.split()[1])
                    literals = map(int, [line.split()[i] for i in range(2,len(line.split())-1)])
                    literals = list(literals) 
                    weight = magic_weight(len(literals),k,'d')
                    if len(literals) == 0:
                        continue
                    if k>=0:
                        result.add_clause(literals,k,'c',weight)         
                    else:
                        result.add_clause([-literals[i] for i in range(len(literals))], len(literals)+k,'c',weight)       
                elif split[1][0] == 'x':
                    literals = []
                    coefs = []
                    split = line.split()
                    for i in range(int((len(split) - 3)/2)):
                        coefs.append(int(split[i*2]))
                        literals.append(int(split[i*2+1][1:]))
                    comparator = split[len(split)-3]
                    k = int(split[len(split)-2])
                    literals, k, coefs, comparator = canonicalize(literals, k, coefs, comparator)
                    weight = magic_weight(len(literals), 1 , 'x') # linear constraints
                    result.add_clause(literals,k,'p',weight,coefs,comparator)
                    #print('add pb constraint. literals= '+repr(literals)+' coefs= '+repr(coefs)+' k= '+repr(k))
                else: # cnf clauses
                    literals = map(int, [line.split()[i] for i in range(0,len(line.split())-1)])
                    literals = list(literals) 
                    coefs = [1 for i in range(len(literals))]
                    weight = magic_weight(len(literals),1,'c')
                    if len(literals) == 0:
                        continue
                    result.add_clause(literals,1,'c',weight, coefs,'>=')
        return result

def magic_weight(n,k,ctype):
    return n
    weight_type = '2'
    negflag = 0
    if k<0:
        negflag = 1
        k = abs(k)
    if weight_type == '1': # n / single inf
        if ctype=='x': return n
        if ctype=='d' or ctype=='g' or ctype=='c':
            return  n * 1.0 * 2**(n-1)/(binom(n-1,abs(k)-1))
    
    if weight_type == '2': # n
        if ctype == 'x': return 2**(n-1)/n
        if ctype == 'c' or ctype == 'd' or ctype == 'g': 
            c = abs(CardinalityFC(n,k)[1])
            print('c='+repr(c))
            return 1/c/n
        return 2**n

    if weight_type == '3': # log(1/prob)
        if ctype=='x': return 1.0
        if ctype=='d' or ctype=='g' or ctype=='c':
            arr = [(binom(n,i)) for i in range(k,n+1)]
            prob = F(sum(arr),(2 ** (n))) * 1.0
            if negflag == 1: prob = 1-prob
            return  -1.0 / math.log2(prob)

    if weight_type == '4': # log(1/prob) / single inf
        if ctype=='x': return 1.0
        if ctype=='d' or ctype=='g' or ctype=='c':
            arr = [(binom(n,i)) for i in range(k,n+1)]
            prob = F(sum(arr),(2 ** (n))) * 1.0
            if prob<1e-5: prob = 1e-5
            if prob>1-1e-5: prob = 1-1e-5
            if negflag == 1: prob = 1-prob
            print('prob= '+repr(prob))
            return  -1.0 / math.log2(prob) * 2**(n-1)/(binom(n-1,abs(k)-1))

    if weight_type == '5': # n / log[1/prob] / single inf
        if ctype=='x': return 1.0
        if ctype=='d' or ctype=='g' or ctype=='c':
            arr = [(binom(n,i)) for i in range(k,n+1)]
            prob = F(sum(arr),(2 ** (n))) * 1.0
            if negflag == 1: prob = 1-prob
            return  -1.0 * n / math.log2(prob) * 2**(n-1)/(binom(n-1,abs(k)-1))

    if weight_type == '6': # n / log[1/prob] 
        if ctype=='x': return 1.0
        if ctype=='d' or ctype=='g' or ctype=='c':
            arr = [(binom(n,i)) for i in range(k,n+1)]
            prob = F(sum(arr),(2 ** (n))) * 1.0
            if negflag == 1: prob = 1-prob
            return  -1.0 * n / math.log2(prob)

    if weight_type == '7': # 1 / (single inf)^2
        if ctype=='x': return 1
        if ctype=='d' or ctype=='g' or ctype=='c':
            return  1.0 * (2**(n-1)/(binom(n-1,abs(k)-1)))**2
    return n

        
 

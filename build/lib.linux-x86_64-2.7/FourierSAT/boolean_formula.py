from collections import defaultdict
from FourierCoefficient import CardinalityFC
import math
import scipy
class Formula:
    def __init__(self):
        self._variables = []
        self._clauses = []
        self._klist = []
        self._weight = []
        self._ctype = []
        self._optsol = 0
    def add_clause(self, literals,k, ctype, weight):
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
                if len(line) == 0 or line[0] == 'c':
                    continue
                if line[0] == 'p':
                    num_vars = int(line.split()[2])
                    result._variables = [i+1 for i in range(num_vars)]
                elif line[0] == 'w':
                    args = line.split()
                    if float(args[2]) == -1:
                        vars[int(args[1])] = result.fresh_variable(1, 1)
                    else:
                        prob = float(args[2])
                        vars[int(args[1])] = result.fresh_variable(1 - prob, prob)
                elif line[0] == 'g':
                    args = line.split()
                    k = int(args[1])
                    weight = magic_weight(num_vars)
                    result._optsol = abs(k)
                    if k>=0:
                        result.add_clause(range(1,num_vars+1),k,'c',weight)#num_vars)#2/(1-c))
                    else:
                        result.add_clause(range(-num_vars,0),num_vars+k,'c', weight)#num_vars)#2/(1+c))
                elif line[0] == 'x':
                    literals = map(int, [line.split()[i] for i in range(1,len(line.split())-1)])
                    literals = list(literals)
                    if len(literals) == 0:
                        continue
                    result.add_clause(literals,0,'x',magic_weight(len(literals)))
                elif line[0] == 'n':
                    literals = map(int, [line.split()[i] for i in range(1,len(line.split())-1)])
                    literals = list(literals)
                    if len(literals) == 0:
                        continue
                    result.add_clause(literals,0,'n',magic_weight(len(literals)))
                elif line[0] == 'd': # cardinality constraints
                    args = line.split()
                    k = int(args[1])
                    literals = map(int, [line.split()[i] for i in range(2,len(line.split())-1)])
                    literals = list(literals) 
                    weight = magic_weight(len(literals))
                    if len(literals) == 0:
                        continue
                    result.add_clause(literals,k,'c',weight)                
                else: # cnf clauses
                    literals = map(int, [line.split()[i] for i in range(0,len(line.split())-1)])
                    literals = list(literals) 
                    weight = magic_weight(len(literals))
                    if len(literals) == 0:
                        continue
                    result.add_clause(literals,1,'c',weight)
        return result

def magic_weight(n):
    return n

        

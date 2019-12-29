FourierSAT is a versatile solver for hybrid Boolean constraints. 

The types of constraints include: CNF (-or), XOR, cardinality constraints and NAE (not all equal).

----------------------------------------------------------------------------------------------------------------------
*Required environment: python with Scipy

*Basic usage:

python usage:
python FourierSAT/FourierSAT.py [DIMACS filepath] --options

command line usage:
python setup.py install
FourierSAT [DIMACS filepath] --options

*Optional parameters:
--timelimit: the time limit (in seconds) for this run. Default: 60
--tolerance: the number of clauses that a solution can violate. Default: 0
--cpus: the number of CPU cores available (the more, the better). Default: 8
-verbose: set it to 1 to output more information

For example:
python FourierSAT.py sample.cnf --timelimit 10 --tolerance 1 --cpus 2 --verbose 1

*Extended DIMACS Format
FourierSAT accepts an extended DIMACS format which can handle CNF, XOR, cardinality constraints and Not-all-equal clauses.

CNF: "[literals] 0"
     eg: clause x_1 or \neg x_2: "1 -2 0"
XOR: "x [literals] 0"
     eg: clause x_1 xor \neg x_2: "x 1 -2 0"
Cardinality constraints: "d [k] [literals] 0"
      k>0 means greater or equal
      k<0 means less or equal
      eg: x_1 + x_2 + x_4 + \neg x_5 >=2: "d 2 1 2 4 -5 0"
      if you want to include a global cardinality constraint (a constriant containing all variables and all the literals are positive), use a line "g [k]". (Note: no '0' at the end of the line!)
      eg: x_1+x_2+x_3+x_4+x_5 <= 2: "g -2"
Not all equal: "n [literals] 0"
      eg: NAE(x_1,x_2,\neg x_3): "n 1 2 -3 0"

*Output: 
-s "solved"/"not-solved in timelimit seconds"+[minimum number of violated clauses]   
-v [solutions]/[the assignment with minimum number of violated clauses found]     

* FourierSAT is a continuous-optimization-based SAT/MaxSAT solver for hybrid Boolean constraints (CNF, XOR, cardinality constraints and not-all-equal) written in Python. 

	- FourierSAT: A Fourier Expansion-Based Algebraic Framework for Solving Hybrid Boolean Constraints
	https://arxiv.org/abs/1912.01032
	(AAAI-2020)

	- A comprehensive journal version: 
	Solving hybrid Boolean constraints in continuous space via multilinear Fourier expansions
	https://akyrillidis.github.io/pubs/Journals/fourierSAT.pdf 
	(AI Journal-2021)

* To achieve better performance, see GradSAT (https://github.com/vardigroup/GradSAT), an extension of FourierSAT. GradSAT is implemented in C++, where BDDs are used for accelerating gradient computations.

	- On Continuous Local BDD-Based Search for Hybrid SAT Solving
	https://arxiv.org/abs/2012.07983
	(AAAI-2021)

* Also check our new work on equiping FourierSAT with GPUs and matrix multiplications with up to 100x accelerated gradient computation: 

	- Massively Parallel Continuous Local Search for Hybrid SAT Solving on GPUs
	https://arxiv.org/abs/2308.15020
        https://github.com/seeder-research/FastFourierSAT

If you have questions or thoughts regarding the tool or this work, please contact zhiwei@rice.edu.

----------------------------------------------------------------------------------------------------------------------
Required environment
-----------------------------------------
 python with Scipy

Basic usage
---------------
python usage:

	python FourierSAT/FourierSAT.py [DIMACS filepath] --options

*Optional parameters:

--timelimit: the time limit (in seconds) for this run. Default: 60

--tolerance: the number of clauses that a solution can violate. Default: 0

--cpus: the number of CPU cores available (the more, the better). Default: 8

--verbose: set it to 1 to output more information

For example:

	FourierSAT sample.cnf --timelimit 10 --tolerance 1 --cpus 2 --verbose 1

Input: Extended DIMACS Format
-------------------------
FourierSAT accepts an extended DIMACS format which can handle CNF, XOR, cardinality constraints and Not-all-equal clauses. MaxSAT instances (.wcnf) and cardinality constraints encoded in pseudo-Boolean format (.opb) are also accepted.

CNF: "[literals] 0"

	eg: clause x_1 or \neg x_2: "1 -2 0"
     
XOR: "x [literals] 0"

     eg: clause x_1 xor \neg x_2: "x 1 -2 0"
     
Cardinality constraints: "d [k] [literals] 0"
      k>0 means greater or equal
      k<0 means less or equal
      
    eg: x_1 + x_2 + x_4 + \neg x_5 >=2: "d 2 1 2 4 -5 0"
    
  Alternatively, you can use the pseudo-Boolean encoding:
  	
	"1 x1 + 1 x2 + 1 x4 - 1 x5 >= 1"
   
  if you want to include a global cardinality constraint (a constriant containing all variables and all the literals are positive), use a line "g [k]". (Note: no '0' at the end of the line!)
  
      eg: x_1+x_2+x_3+x_4+x_5 <= 2: "g -2"
      
Not all equal: "n [literals] 0"

      eg: NAE(x_1,x_2,\neg x_3): "n 1 2 -3 0"

Output
-------
	-s "solved"/"not-solved in timelimit seconds"+[minimum number of violated clauses]   
	-v [solutions]/[the assignment with minimum number of violated clauses found]    
	-o [the cost of the best solution found so far (the number of violated constraints)] (Only for MaxSAT mode)

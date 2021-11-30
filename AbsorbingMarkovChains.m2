
newPackage(
    "AbsorbingMarkovChains",
    Version => "0.1",
    Date => "Fall 2021",
    Headline => "Finding the fundamental matrix and time to absorption for an Absorbing Markov Chain",
    Authors => {{ Name => "Geeth M.", Email => "gmahagamage1206@lions.piedmont.edu"}}
    )

export {"fundamentalMatrix","timeToAbsorption","absorptionProbabilities"}

-* Code section *-

TransitionMatrix = new SelfInitializingType of Matrix

transitionMatrix := method()
transitionMatrix Matrix := A -> (
    if numRows A!= numColumns A then
       error "a square matrix is expected";
    if not all(sum entries transpose A, x -> x==1) then
       error "The rows don't add upto 1";
    TransitionMatrix A )


getQandR = method()

getQandR TransitionMatrix := A -> (
    div := position(entries A, row -> any(row, x -> x == 1));
    n := numRows A;
    Q := A_{0..(div - 1)}^{0..(div - 1)};
    R := A_{div..(n - 1)}^{0..(div - 1)};    
    (Q, R))

getQandR Matrix := A -> getQandR transitionMatrix A


fundamentalMatrix = method()

fundamentalMatrix TransitionMatrix := A ->(
    (Q,R) := getQandR A;
    I := id_(RR^(numRows Q));
    inverse (I-Q))

fundamentalMatrix Matrix := A -> fundamentalMatrix transitionMatrix A




timeToAbsorption = method()

timeToAbsorption TransitionMatrix := A -> (
    (Q,R) := getQandR A;
    N := fundamentalMatrix A;
    c := transpose matrix { for i from 1 to (numRows Q) list 1};
    N*c)

timeToAbsorption Matrix := A -> timeToAbsorption transitionMatrix A


absorptionProbabilities = method()

absorptionProbabilities TransitionMatrix := A -> (
    (Q,R) := getQandR A;
    N := fundamentalMatrix A ;
    N*R)


absorptionProbabilities Matrix := A -> absorptionProbabilities transitionMatrix A

beginDocumentation()


doc ///
  Key 
    AbsorbingMarkovChains
  Headline 
  --headline
  Description
   Text
    --- a description of the package--

///

doc ///
  Key
    getQandR
    (getQandR, Matrix)
  Headline
    get Q and R Matrices  
  Usage
    getQandR(A)
  Description
   Text
    This is an intermediary function which will take in a nxn matrix
    in its canonical form and will break it into Q and R sections 
    used for the next few functions
   Example
     debug AbsorbingMarkovChains
     A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}}
     getQandR A

   
/// 

doc ///
  Key
    fundamentalMatrix
    (fundamentalMatrix, Matrix)
  Usage
    fundamentalMatrix(A)
  Description
   Text
    This function takes in a nxn matrix in its canonical form and 
    outputs its fundmental matrix 
   Example
     A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}}
     fundamentalMatrix A




/// 

doc ///
  Key
    timeToAbsorption
    (timeToAbsorption, Matrix)
  Usage
    timeToAbsorption(A)
  Description
   Text
    This function takes in a nxn matrix in its canonical form and 
    outputs its Time To Absorption matrix
   Example
     A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}}
     timeToAbsorption A




/// 

doc ///
  Key
    absorptionProbabilities
    (absorptionProbabilities, Matrix)
  Usage
    absorptionProbabilities(A)
  Description
   Text
     This function takes in a nxn matrix in its canonical form and 
    outputs its Absorption Probabilities matrix
   Example
     A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}}
     absorptionProbabilities A




/// 


TEST ///
A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}} 
debug AbsorbingMarkovChains
assert(getQandR A == (matrix{{0,1/2,0},{1/2,0,1/2},{0,1/2,0}},matrix{{1/2,0},{0,0},{0,1/2}}))
///


TEST ///
A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}} 
assert zero clean(1e-15, (fundamentalMatrix A - matrix{{1.5,1,0.5},{1,2,1},{0.5,1,1.5}}))
///

TEST ///
A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}} 
assert zero clean(1e-15, (timeToAbsorption A - matrix{{3},{4},{3}}))
///

TEST ///
A = matrix{{0,1/2,0,1/2,0},{1/2,0,1/2,0,0},{0,1/2,0,0,1/2},{0,0,0,1,0},{0,0,0,0,1}} 
assert zero clean(1e-15, (absorptionProbabilities A - matrix{{0.75,0.25},{0.5,0.5},{0.25,0.75}}))
///
 
















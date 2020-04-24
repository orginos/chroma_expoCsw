Tests for 6x6 Traceles Hermitian Exponential
---

24/4/2020
- Here are 2 simple files to test the recursive algorithm to compute the exponential of the matrices for the clover term. 

- There is a _python_ implementation which contains both the recursive method and a direct calculation, using whatever algorithm implemented in _numpy_, to compare against.

- The _C++_ implementation uses the same data layout for the matrix as the one used in _Chroma_. The results match the ones from the python code.

- The _C++_ has been added to the ```lib/actions/ferm/linop/expo_clover_term_qdp.h``` file, using the correct datatypes. Please check if some things could have been done more efficiently. 

- The input matrix still has to be fixed
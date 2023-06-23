# Topological elements
It is possible to compare datasets based on C3 and S(k) counts. 

In our paper, we introduced the concept of a 3-vertex-cycle, or 'C3', which is a cycle comprising three links that all share at least one common tissue. More formally, a triplet of nodes `(A, B, D)` forms a C3 if these nodes share a common tissue across their interactions `(A, B) ∈ Eₜ`, `(B, D) ∈ Eₜ`, and `(A, D) ∈ Eₜ`. We say that a C3 `(A, B, D)` includes a link `e` if `e` is one of its three links. `Eₜ` is Hi-C interaction graph of one tissue `t`. 

The 'degree' of a C3, denoted as `DegC3(A, B, D)`, is defined as the number of tissues that share the same interactions that form the C3. If the degree of a triplet of nodes is more than 0 (`DegC3(A, B, D) > 0`), these nodes form a C3.

A link that is part of at least `k` different C3 cycles is defined as a 'support', denoted as `S(k)`. The 'degree' of a support link, symbolized as `DegSupport(A, B)`, is the count of distinct C3 cycles that include the link, that is, `DegSupport(A, B) = |{ D | D ∈ V and (A, B, D) is C3 }|`.

## C3 and S(k) calculation 
Notebook compareDatasets.ipynb demonstrates the calculation of C3 and S(k), and creates images to compare C3 and S(k) of 3 datasets.
It also creates a csv file with C3 and S(k) counts for 3 datasets
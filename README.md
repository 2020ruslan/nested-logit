

![Alt text](ONL.png?raw=true "Optimal Nested Logit")
# Introduction
This package implements the nested logit strucutre learning algirthim as disscussed in Aboutaleb et al. 2020.

This work is about developing an estimation procedure for nested logit models that optimizes over
the nesting structure in addition to the model parameters. Current estimation practices require
an a priori specification of a nesting structure. We formulate the problem of learning an optimal
nesting structure as a mixed integer nonlinear programming (MINLP) optimization problem and
solve it using a variant of the linear outer approximation algorithm. We demonstrate that it is
indeed possible to recover the nesting structure directly from the data by applying our method to
synthetic and real datasets.

For more details: https://dspace.mit.edu/bitstream/handle/1721.1/123208/1129597025-MIT.pdf?sequence=1&isAllowed=y 

# Installation 


```julia
julia> Pkg.add("JuMP")
```

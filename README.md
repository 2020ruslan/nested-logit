

![Alt text](ONL.png?raw=true "Optimal Nested Logit")
# Introduction
This package implements the nested logit strucutre learning algirthim as proposed in Aboutaleb, Ben-Akiva, and Jaillet (2020) https://arxiv.org/abs/2008.08048.

This work is about developing an estimation procedure for nested logit models that optimizes over
the nesting structure in addition to the model parameters. Current estimation practices require
an a priori specification of a nesting structure. We formulate the problem of learning an optimal
nesting structure as a mixed integer nonlinear programming (MINLP) optimization problem and
solve it using a variant of the linear outer approximation algorithm. We demonstrate that it is
indeed possible to recover the nesting structure directly from the data by applying our method to
synthetic and real datasets.


# Installation 

see Installation.pdf
```julia
julia> Pkg.add("JuMP")
```
# Loading data from CSV files 

```julia
using CSVFiles, DataFrames, CSV


main_dataframe=DataFrame(load("Filedirectory/file.csv"))

# filter records by row if needed
main_dataframe=filter(row -> row[:MODE] <= 6, main_dataframe)
main_dataframe=filter(row -> row[:MODE] >0, main_dataframe)

n,~=size(main_dataframe)

#transform columns
main_dataframe[:CAR_TT]=main_dataframe[:CAR_TT]./60
main_dataframe[:PT_TT]=main_dataframe[:PT_TT]./60
#main_dataframe=CSV.read("CTPS.csv")

#observed choices
Y=main_dataframe[:MODE]

#attributes
X=main_dataframe


#split into train and validation
train_validation_ratio=0.75
idx = shuffle(1:n)
train_idx = view(idx, 1:floor(Int, train_validation_ratio*n))
test_idx = view(idx, (floor(Int, train_validation_ratio*n)+1):n)


Y_train=Y[train_idx]
X_train=X[train_idx,:]

Y_validation=Y[validation_idx]
X_validation=X[validation_idx,:]
```


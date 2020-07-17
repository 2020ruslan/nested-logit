
#Import packages
using JuMP, ForwardDiff, Gurobi, Distributions, LinearAlgebra, Random
using LightGraphs,Plots, TikzGraphs
using Ipopt

#if you are using KNITRO you need to import it here: using KNITRO
#you also need to comment line 989 and uncomment line 990 in NLhelperfunctions.jl
#Ipopt should work fine and hey... it is open source

#Include nested logit helper functions:
include("NLhelperfunctions.jl")


#In this example we generate monte carlo data according to a known nesting tree and systematic specficiation and then recover that tree from the data
# We consider a case with four alternatives where alt 1 and alt 2 are in a nest and alt 3 and alt 4 are in another nest.

num_alternatives=4
sample_size=3000

#for 4 alternatives there can be atmost 4-2 = 2 nest nodes excluding the root node. see theorem in thesis. You can specify a tree for the monte carlo
# data generation by specifiying an ancestry vector:
tree_ancestry=[0 1 1 2 2 3 3]
#This is 1 + 2 + 4 vector denoting the parents of the root, the two nests and each of the four alternatives respectively
# The root has no parent so we enter 0, nests 1 and nest 2 have the root as their parent so we enter 2 and finally the parents
#for alt 1 and alt 2 are nest 2 and alt 3 and alt 4 are nest 3 (we count here the root as a nest)
#For real datasets you don't have to worry about understanding or using this notation.

#generate a matrix represenation for this tree
tree_matrix=construct_graph_from_ancestry(tree_ancestry)

#Nest we need to specify the systematic component for this we declare a specficiation matrix S
#S will have dimensions as follows: num_alternatives X num_parameters. S_ij =1 imples that the utility for alternative i includes parameter j multiplied by X_j
#where X is the data matrix.

#In this simple example we consider a specification with only alternative specific constant (ASC)
# we specify our the paramters that enter the systematic specification

true_betas=[0 0.5 0.5 0.5]
num_parameters=length(true_betas)
#here we have four ASCs.

S=zeros(num_alternatives,num_parameters)
for i=1:num_alternatives
    S[i,i]=1
end

#Alterantively you may code the utility functions directly in the helper functions file see line 328 to 346 for an example
#if you decide to hard code the utility you need to comment lines 348 to 350


#Next we need to specify the scale parameters for the nests
#the root has scale one and is fixed by default. The two nests have scales 2 and 4 respectively.
true_scales=[1 2 4]

#Generate the montecarlo data
train_to_validation_ratio=2/3
availability_prob=0.85
A_train,C_train,X_train,Y_train,A_val,C_val,X_val,Y_val=generate_montecarlo_data(sample_size,S,tree_matrix,true_betas,true_scales,availability_prob,train_to_validation_ratio)

#run the nested logit strucutre learning algorithm
max_iterations=5
#this is the maximum number of hyperlanes used to apprpximate the feasibility set. See thesis for details. I recommend starting with 5 and increasing it to 10 with more alternatives or smaller sample size
optimal_trees,optimal_ll,optimal_params,optimal_std_errs,counts=solve_NLSLP(S,X_train,A_train,Y_train,X_val,A_val,Y_val,max_iterations)


#Done!!!

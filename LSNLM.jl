#include("NL_estimation.jl")
using JuMP, Ipopt, KNITRO, ForwardDiff, Gurobi, Distributions, Juniper, Cbc

using LightGraphs
using GraphPlot
using Plots

using QuantEcon: meshgrid
n=1000 #sample size
epsilon=1e-300 #to prevent log(0)
m=8; #number of alternatives
dim1=(1+(m-2));
dim2=(1+(m-2)+m);
index_shift=dim1*dim2
p=2*m;

delta=+0.001;
mu_norm=1;
mu_upper=10;

#######################################################################################################
#Helper functions
feval(fn_str, args...) = eval(parse(fn_str))(args...)

function multiply(u,v,S,X,A,x,b,mu)
    if u==0
        return 0
    else
        return u*feval(v,S,X,A,x,b,mu);
    end
end

####################################################
####################################################33333333
function print_edges(x)
    A=convert(Array{Int64},vcat(round.(reshape(x,dim1,dim2),2),zeros(m,dim2)))
    Tree=DiGraph(A)
    for e in edges(Tree)
       println(e)
      end
end

function construct_ancestry(x)
    nodes=zeros(1+(m-2)+m,1);
    for i=2:1+(m-2)+m
        parent=find(y-> (y==1),x[:,i]);
        if length(parent)>0
            nodes[i]=find(y-> (y==1),x[:,i])[1];
        else
            nodes[i] =0
        end
    end
    nodes=convert(Array{Int64},nodes)
    return nodes
end
function construct_adjacency(nodes)
    x=zeros(dim1,dim2)
    for i=1:dim2
        if nodes[i]>0
            x[nodes[i],i]=1
        end
    end
    return x
end



#Define Utilities
for i=1:m
    str=string("V_",i,"(S,X,b) = vecdot(b.*S[",i,",:],X)");
    eval(Meta.parse(str))
end

#Define Inclusive Values
for i=1:1+(m-2)
    str=string("IV_",(i),"(S,X,A,x,b,mu) = log(epsilon+");

    for j=1+(m-2)+1:(m-2)+m+1
        str=string(str,"A[",(j-((m-2)+1)),"]*","x[",(i),",",(j),"]","*exp(mu[",(i),"]*V_",(j-((m-2)+1)),"(S,X,b))+");
    end
    for j=2:(m-2)+1
        if j==i
            continue
        end
        arg1=string("x[",(i),",",(j),"]");
        arg2=string("\"","IV_",(j),"\"");
        str=string(str,arg1,"*exp(mu[",(i),"]*multiply(",arg1,",", arg2,",S,X,A,x,b,mu))+");
    end
    str=string(str,"0)/(mu[",(i),"]+epsilon)");
    #println(str)
    eval(Meta.parse(str))
end

#Define Choice probabilities
function Probability(nodes,i,j,S,X,A,x,b,mu)
    parent=nodes[j]
    temp=mu[nodes[j]]*feval(string("V_",j-1-(m-2)),S,X[i,:],b)
    anscestor=nodes[parent]
    while anscestor>0
        temp=temp+(mu[anscestor]-mu[parent])*feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)
        parent=anscestor
        anscestor=nodes[parent]
    end
    temp=temp-mu[parent]*feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)
    return temp
end

function Probability_unnormalized(nodes,i,j,S,X,A,x,b,mu)
    parent=nodes[j]
    temp=mu[nodes[j]]*feval(string("V_",j-1-(m-2)),S,X[i,:],b)
    anscestor=nodes[parent]
    while anscestor>0
        temp=temp+(mu[anscestor]-mu[parent])*feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)
        parent=anscestor
        anscestor=nodes[parent]
    end
    return temp
end

function construct_probability_matrix(S,X,A,x,b,mu)
    nodes=construct_ancestry(x);
    Q=zeros(n,m);
    for i=1:n
        for j=1+(m-2)+1:1+(m-2)+m
            Q[i,j-1-(m-2)]=A[i,j-1-(m-2)]*exp(Probability(nodes,i,j,S,X,A,x,b,mu))
        end
    end
    return Q
end

function construct_choice_matrix(Q)
    Y=m-sum(rand(n,1).<cumsum(Q,2),2)+1;
    C=zeros(n,m);
    for i=1:n
        for j=1:m
            if j==Y[i,1]
                C[i,j]=1;
                break
            end
        end
    end
    return C
end

function likelihood(b,mu)
    likelihood=0
    nodes=construct_ancestry(x_t);
    Q=zeros(n,m);
    @inbounds @simd for i=1:n
        for j=1+(m-2)+1:1+(m-2)+m
            #println(i,n)
            #println("beta",b)
            #println()
            likelihood+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(Probability(nodes,i,j,S_t,X,A,x_t,b,mu))
        end
    end
    return likelihood
end



function path_from_root(j,nodes)
    path_j=[]
    node=j
    push!(path_j,node)
    while node>0
        node=nodes[node]
        push!(path_j,node)
        if length(path_j)>=m
            break
        end
    end
    return reverse(path_j[1:end-1])
end

function path_has_edge(path,i,j)
    for k=1:length(path-1)
        if i==path[k] && j==path[k+1]
            return true
        end
    end
    return false
end

function path_has_node(path,i)
    for k=1:length(path)
        if i==path[k]
            return true,k
        end
    end
    return false,-1
end



#root r to alternative j
function gradient_ra(j,b,mu,x,S,A,X)
    grad=0;
    @simd for i=1:n
        util=feval(string("V_",j),S,X[i,:],b)
        utilr=feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
        grad+=C[i,j]*A[i,j]*mu[1]*(util-utilr)
        for k=1:m
            grad+=C[i,k]*A[i,k]*(-exp(mu[1]*(util-utilr)))
        end
    end
    return grad
end

function gradient_ra_ra(j1,j2,b,mu,x,S,A,X)
    grad=0;
    @simd for i=1:n
        util1=feval(string("V_",j1),S,X[i,:],b)
        util2=feval(string("V_",j2),S,X[i,:],b)
        utilr=feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
        for k=1:m
            grad+=C[i,k]*A[i,k]*(exp(mu[1]*(util1+util2-2*utilr)))/mu[1]
        end
    end
    return grad
end



function gradient_rb_rb(k1,k2,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x)
    @simd for j=1+(m-2)+1:1+(m-2)+m
        for i=1:n
            norm_const=feval(string("IV_",k1),S,X[i,:],A[i,:],x,b,mu)+feval(string("IV_",k2),S,X[i,:],A[i,:],x,b,mu)-2*feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
            grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(exp(mu[1]*norm_const)/mu[1])
        end
    end
    return grad
end

function gradient_rb_ra(k,j1,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x)
    @simd for j=1+(m-2)+1:1+(m-2)+m
        for i=1:n
            norm_const=feval(string("V_",j1),S,X[i,:],b)+feval(string("IV_",k),S,X[i,:],A[i,:],x,b,mu)-2*feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
            grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(exp(mu[1]*norm_const)/mu[1])
        end
    end
    return grad
end


function gradient_bb(k,l,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x)
    @simd for j=1+(m-2)+1:1+(m-2)+m
        path=path_from_root(j,nodes)
        if path_has_edge(path,k,l)
            for i=1:n
                parent=nodes[l]
                ancestor=nodes[parent]
                norm_const=((mu[ancestor]-mu[parent])/mu[parent])*exp(mu[parent]*(feval(string("IV_",l),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)))
                grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(Probability_unnormalized(nodes,i,j,S_t,X,A,x_t,b,mu)+norm_const)
            end
        elseif path_has_node(path,k)
            for i=1:n
                parent=k
                ancestor=nodes[parent]
                norm_const=((mu[ancestor]-mu[parent])/mu[parent])*exp(mu[parent]*(feval(string("IV_",l),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)))
                grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(norm_const)
            end

        else
            grad+=0
        end

    end
    return grad
end

function gradient_bb_bb(k1,l1,k2,l2,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x)
    if k1==k2 && l2==l2
        for j=1+(m-2)+1:1+(m-2)+m
            path=path_from_root(j,nodes)
            if path_has_edge(path,k1,l1)&&path_has_edge(path,k2,l2)
                #println("path from root to ", j ," through ", k)
                for i=1:n
                    parent1=nodes[l1]
                    ancestor1=nodes[parent1]

                    util1=feval(string("IV_",l1),S,X[i,:],A[i,:],x,b,mu)
                    util2=feval(string("IV_",parent1),S,X[i,:],A[i,:],x,b,mu)

                    norm_const1=((mu[ancestor1]-mu[parent1])/mu[parent1])*exp(mu[parent1]*(util1-util2))
                    norm_const2=-1*((mu[ancestor1]-mu[parent1])/mu[parent1])*exp(mu[parent1]*(2*util1-2*util2))
                    grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(norm_const1+norm_const2)

                end
            end
        end

    else
        for j=1+(m-2)+1:1+(m-2)+m
            path=path_from_root(j,nodes)
            if path_has_edge(path,k1,l1)&&path_has_edge(path,k2,l2)
                #println("path from root to ", j ," through ", k)
                for i=1:n
                    parent1=nodes[l1]
                    ancestor1=nodes[parent1]

                    parent2=nodes[l2]
                    ancestor2=nodes[parent2]

                    norm_const1=((mu[ancestor1]-mu[parent1])/mu[parent1])*exp(mu[parent1]*(feval(string("IV_",l1),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",parent1),S,X[i,:],A[i,:],x,b,mu)))
                    norm_const2=((mu[ancestor2]-mu[parent2])/mu[parent2])*exp(mu[parent2]*(feval(string("IV_",l2),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",parent2),S,X[i,:],A[i,:],x,b,mu)))
                    grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(Probability_unnormalized(nodes,i,j,S_t,X,A,x_t,b,mu)+norm_const1+norm_const2)
                end
            end
        end
    end

    return grad
end

function gradient_bb_rb(k1,l1,k2,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x)
    @simd for j=1+(m-2)+1:1+(m-2)+m
        path=path_from_root(j,nodes)
        if path_has_edge(path,k1,l1)&&path_has_edge(path,1,k2)
            #println("path from root to ", j ," through ", k)
            for i=1:n
                parent=nodes[l1]
                ancestor=nodes[parent]
                norm_const=((mu[ancestor]-mu[parent])/mu[parent])*exp(mu[parent]*(feval(string("IV_",l1),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)))
                grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(Probability_unnormalized(nodes,i,j,S_t,X,A,x_t,b,mu)+norm_const)
            end
        end
    end
    return grad
end


function gradient_ba_rb(k1,j1,k2,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x)
    path=path_from_root(j1+1+(m-2),nodes)
    if path_has_edge(path,k1,j1+1+(m-2))&&path_has_edge(path,1,k2)
        for i=1:n
            parent=nodes[j1+1+(m-2)]
            ancestor=nodes[parent]
            norm_const=((mu[ancestor]-mu[parent])/mu[parent])*exp(mu[parent]*(feval(string("V_",j1),S,X[i,:],b)-feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)))
            grad+=C[i,j1]*A[i,j1]*(Probability_unnormalized(nodes,i,j1+1+(m-2),S_t,X,A,x_t,b,mu)+norm_const)
        end
    end
    return grad
end





function gradient_ba_ba(k1,j1,k2,j2,b,mu,x,S,A,X)
    grad=0
    if k1==k2 && j1==j2
        nodes=construct_ancestry(x)
        path=path_from_root(j1+1+(m-2),nodes)
        if path_has_edge(path,k1,j1+1+(m-2))
            for i=1:n
                parent=nodes[j1+1+(m-2)]
                ancestor=nodes[parent]

                util1=feval(string("V_",j1),S,X[i,:],b)
                util2=feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)

                norm_const1=((mu[ancestor]-mu[parent])/mu[parent])*exp(mu[parent]*(util1-util2))
                norm_const2=-1*((mu[ancestor]-mu[parent])/mu[parent])*exp(mu[parent]*(2*util1-2*util2))
                grad+=C[i,j1]*A[i,j1]*(norm_const1+norm_const2)
            end
        end
    end
    return grad

end

function gradient_bb_ba(k1,l1,k2,j1,b,mu,x,S,A,X)
    grad=0
    nodes=construct_ancestry(x);
    path=path_from_root(j1+1+(m-2),nodes)
    if path_has_edge(path,k2,j1+1+(m-2))&&path_has_edge(path,k1,l1)
        #println("path from root to ", j ," through ", k)
        for i=1:n
            parent1=nodes[l1]
            ancestor1=nodes[parent1]

            parent2=nodes[j1+1+(m-2)]
            ancestor2=nodes[parent2]

            norm_const1=((mu[ancestor1]-mu[parent1])/mu[parent1])*exp(mu[parent1]*(feval(string("IV_",l1),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",parent1),S,X[i,:],A[i,:],x,b,mu)))
            norm_const2=((mu[ancestor2]-mu[parent2])/mu[parent2])*exp(mu[parent2]*(feval(string("V_",j1),S,X[i,:],b)-feval(string("IV_",parent2),S,X[i,:],A[i,:],x,b,mu)))
            grad+=C[i,j1]*A[i,j1]*(Probability_unnormalized(nodes,i,j1+1+(m-2),S_t,X,A,x_t,b,mu)+norm_const1+norm_const2)
        end
    end
    return grad
end
#root r to nest k

#nest k to nest l



function define_gradients()
    dim1=(1+(m-2));
    dim2=(1+(m-2)+m);
    grad=zeros(dim1,dim2);

    for i=1:1
        for j=2:dim1
            func_name=string("gradient_rb_",j)
            str=string(func_name,"(y) =gradient_rb(",j,",y[1:p],y[p+1:end],x_t,S_t,A,X)")
            eval(Meta.parse(str))
            str=string("fd_",func_name,"= y->ForwardDiff.gradient(",func_name,",y)")
            eval(Meta.parse(str))
        end
        for j=dim1+1:dim2
            func_name=string("gradient_ra_",j-1-(m-2))
            str=string(func_name,"(y) =gradient_ra(",j-1-(m-2),",y[1:p],y[p+1:end],x_t,S_t,A,X)")
            eval(Meta.parse(str))
            str=string("fd_",func_name,"= y->ForwardDiff.gradient(",func_name,",y)")
            eval(Meta.parse(str))
        end
    end

    for j=2:dim1
        for k=2:dim1
            if j==k
                continue
            else
                func_name=string("gradient_bb_",j,"_",k)
                str=string(func_name,"(y) =gradient_bb(",j,",",k,",y[1:p],y[p+1:end],x_t,S_t,A,X)")
                eval(Meta.parse(str))
                str=string("fd_",func_name,"= y->ForwardDiff.gradient(",func_name,",y)")
                eval(Meta.parse(str))
            end
        end
    end

    for i=2:dim1
        for j=1:m
            func_name=string("gradient_ba_",i,"_",j)
            str=string(func_name,"(y) =gradient_ba(",i,",",j,",y[1:p],y[p+1:end],x_t,S_t,A,X)")
            eval(Meta.parse(str))
            str=string("fd_",func_name,"= y->ForwardDiff.gradient(",func_name,",y)")
            eval(Meta.parse(str))
        end
    end

end

define_gradients()

function generalized_utility(node,person,β,μ,x,S,A,X)
    if node<=dim1
        #return feval(string("IV_",node),S,X[person,:],A[person,:],x,β,μ)
        return max(feval(string("IV_",node),S,X[person,:],A[person,:],x,β,μ),0)
    else
        return feval(string("V_",node-dim1),S,X[person,:],β)
    end
end

#δΓa wrt x_bc
# a: nest
# b: nest
# c: node
function δΓ(a,b,c,person,β,μ,x,S,A,X)
    if a==b
        util_diff=generalized_utility(c,person,β,μ,x,S,A,X)-generalized_utility(a,person,β,μ,x,S,A,X)
        return (1/μ[a])*exp(μ[a]*util_diff)
    else
        grad=0
        for k=2:dim1
            if x[a,k]>0
                util_diff=generalized_utility(k,person,β,μ,x,S,A,X)-generalized_utility(a,person,β,μ,x,S,A,X)
                grad+=x[a,k]*exp(μ[a]*util_diff)*δΓ(k,b,c,person,β,μ,x,S,A,X)
                #println(exp(μ[a]*util_diff))
                #println("k ", k)
            end
        end
        return grad
    end
end

function accumulate_δΓ(nodes,derivative_nest,path_alternative,derivative_alternative,person,S,X,A,x,b,mu)
    temp=0
    parent=nodes[path_alternative]
    anscestor=nodes[parent]
    while anscestor>0
        temp=temp+(mu[anscestor]-mu[parent])*δΓ(parent,derivative_nest,derivative_alternative,person,b,mu,x,S,A,X)
        parent=anscestor
        anscestor=nodes[parent]
    end
    temp=temp-mu[1]*δΓ(1,derivative_nest,derivative_alternative,person,b,mu,x,S,A,X)
    return temp
end


function gradient_ba(nest,alternative,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x);
    for i=1:n
        for j=dim1+1:dim1+m
            if j==alternative
                if C[i,j-dim1]*A[i,j-dim1]>0
                    grad+=C[i,j-dim1]*A[i,j-dim1]*accumulate_δΓ(nodes,nest,alternative,alternative,i,S,X,A,x,b,mu)
                end
                #parent=nodes[nest]
                flag,~=path_has_node(path_from_root(j,nodes),nest)
                if flag
                    utilr=feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
                    adjusted_ancestry=nodes
                    adjusted_ancestry[alternative]=nest
                    grad+=C[i,j-dim1]*A[i,j-dim1]*(Probability_unnormalized(adjusted_ancestry,i,j,S_t,X,A,x,b,mu)-mu[1]*utilr)
                end
            else
                if C[i,j-dim1]*A[i,j-dim1]>0
                    grad+=C[i,j-dim1]*A[i,j-dim1]*accumulate_δΓ(nodes,nest,alternative,alternative,i,S,X,A,x,b,mu)
                end
            end
        end
    end
    return grad
end

function gradient_bb(nest_origin,nest_destination,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x);
    for i=1:n
        for j=dim1+1:dim1+m
            #println(i,j)
            if C[i,j-dim1]*A[i,j-dim1]>0
                grad+=C[i,j-dim1]*A[i,j-dim1]*accumulate_δΓ(nodes,nest_origin,j,nest_destination,i,S,X,A,x,b,mu)
            end

            path_from_root_to_alternative=path_from_root(j,nodes)
            flag1,position1=path_has_node(path_from_root_to_alternative,nest_origin)
            flag2,position2=path_has_node(path_from_root_to_alternative,nest_destination)
            if flag1 && flag2 && position1<position2
                utilr=feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
                adjusted_ancestry=nodes
                adjusted_ancestry[nest_destination]=nest_origin
                grad+=C[i,j-dim1]*A[i,j-dim1]*(Probability_unnormalized(adjusted_ancestry,i,j,S_t,X,A,x_t,b,mu)-mu[1]*utilr)
            end
        end
    end
    return grad
end


function gradient_rb(nest,b,mu,x,S,A,X)
    grad=0;
    nodes=construct_ancestry(x)
    @simd for i=1:n
        utilk=feval(string("IV_",nest),S,X[i,:],A[i,:],x,b,mu)
        utilr=feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
        norm_const=utilk-utilr
        for j=1+(m-2)+1:1+(m-2)+m
            flag,~=path_has_node(path_from_root(j,nodes),nest)
            if flag
                adjusted_ancestry=nodes
                adjusted_ancestry[nest]=0
                C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(Probability_unnormalized(adjusted_ancestry,i,j,S_t,X,A,x_t,b,mu)-mu[1]utilr-exp(mu[1]*norm_const))
            else
                grad+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(-1*exp(mu[1]*norm_const))
            end
        end
    end
    return grad
end

function gradient_x(b,mu,x,S,A,X)
    dim1=(1+(m-2));
    dim2=(1+(m-2)+m);
    grad=zeros(dim1,dim2);

    for i=1:1
        for j=2:dim1
            grad[i,j]=gradient_rb(j,b,mu,x,S,A,X)
        end
        for j=dim1+1:dim2
            grad[i,j]=gradient_ra(j-1-(m-2),b,mu,x,S,A,X)
        end
    end

    for j=2:dim1
        for k=2:dim1
            if j==k
                continue
            else
                grad[j,k]=gradient_bb(j,k,b,mu,x,S,A,X)
            end
        end
    end

    for i=2:dim1
        for j=1+dim1:m+dim1
            grad[i,j]=gradient_ba(i,j,b,mu,x,S,A,X)
        end
    end

    return -1*grad
end
#δΓ_a1/δ_(a2,j)
#function δΓ(a1,a2,j,i,b,mu,x,S,A,X)
#    if a1==1
#        grad=0
#        for k=2:dim1
#            norm_const=feval(string("IV_",k),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",1),S,X[i,:],A[i,:],x,b,mu)
#            grad+=x[1,k]*exp(mu[1]*norm_const)*δΓ(k,a2,j,i,b,mu,x,S,A,X)
#        end
#        return grad
#    end
#    if a1==a2
#        norm_const=feval(string("V_",j),S,X[i,:],b)-feval(string("IV_",a2),S,X[i,:],A[i,:],x,b,mu)
#        return (1/mu[a2])*exp(mu[a2]*(norm_const))
#    else
#        norm_const=feval(string("IV_",a2),S,X[i,:],A[i,:],x,b,mu)-feval(string("IV_",a1),S,X[i,:],A[i,:],x,b,mu)
#        return x[a1,a2]*exp(mu[a1]*(norm_const))*δΓ(a2,a2,j,i,b,mu,x,S,A,X)
#    end
#end



function Probability(nodes,i,j,S,X,A,x,b,mu)
    parent=nodes[j]
    temp=mu[nodes[j]]*feval(string("V_",j-1-(m-2)),S,X[i,:],b)
    anscestor=nodes[parent]
    while anscestor>0
        temp=temp+(mu[anscestor]-mu[parent])*feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)
        parent=anscestor
        anscestor=nodes[parent]
    end
    temp=temp-mu[parent]*feval(string("IV_",parent),S,X[i,:],A[i,:],x,b,mu)
    return temp
end


function likelihood_x(b,mu,x)
    likelihood=0
    nodes=construct_ancestry(x);
    Q=zeros(n,m);
    @simd for i=1:n
        for j=1+(m-2)+1:1+(m-2)+m
            likelihood+=C[i,j-1-(m-2)]*A[i,j-1-(m-2)]*(Probability(nodes,i,j,S_t,X,A,x,b,mu))
        end
    end
    return likelihood
end

function likelihood_2(y)
    return -1*likelihood(y[1:2*m],y[(2*m+1):end])
end

function f_x(y)
    dim1=(1+(m-2));
    dim2=(1+(m-2)+m);
    index_shift=dim1*dim2;
	return -1*likelihood_x(y[1:2*m],y[(2*m+1):(end-index_shift)],reshape(y[(end-index_shift+1):end],dim1,dim2))
end

function f(y)
	return -1*likelihood(y[1:2*m],y[(2*m+1):(end)])
end

function fouter_JuMP(y...)
	return -1*likelihood(y[1:2*m],y[(2*m+1):end])
end

function fouter_x_JuMP(y...)
    dim1=(1+(m-2));
    dim2=(1+(m-2)+m);
    index_shift=dim1*dim2;
	return -1*likelihood_x(y[1:2*m],y[(2*m+1):(end-index_shift)],convert(Array{Int64},reshape(y[(end-index_shift+1):end],dim1,dim2)))
end

function solve_NL_subproblem()
    println("Solving NL subproblem")
    #model=Model(solver=CplexSolver(CPXPARAM_OptimalityTarget=3))
    #model= Model(solver=KnitroSolver(algorithm=4,maxit=80,par_numthreads=3, convex=0))
    #ipopt = IpoptSolver(print_level=0); cbc = CbcSolver()
    #model=  Model(solver=JuniperSolver(ipopt, mip_solver=cbc))
    model= Model(solver=KnitroSolver(algorithm=2,maxit=5000))
    #model=Model()
    #model = Model(solver=IpoptSolver(max_iter=100))
    JuMP.register(model, :fouter_JuMP, 2*m+1+(m-2),fouter_JuMP, autodiff=true)
    @variable(model, y[1:2*m+1+(m-2)])

    #for i=1
        #setvalue(y[i], warm_start_vector[i])
    #end


    for i=1:1+(m-2)
        for j=1:1+(m-2)
            if x_t[i,j]==1
                #println(i,j)
                @constraint(model,y[2*m+j]>=delta+y[2*m+i])
                #@constraint(model,y[2*m+j]-1==y[2*m+i])
            end
        end
    end

    #add_nest_constraints(x_t)

    for j=2:1+(m-2)
        @constraint(model,y[2*m+j]>=mu_norm)
        @constraint(model,y[2*m+j]<=mu_upper)
    end

    for j=1:m
        @constraint(model,y[j]<=10)
    end

    for j=m+1:2*m
        @constraint(model,y[j]<=1)
    end


    @constraint(model,y[2*m+1]==mu_norm)
    @constraint(model,y[m]==0)

    #@NLobjective(model,:Max,fouter)
    objexpr = Expr(:call, :fouter_JuMP, y...)
    JuMP.setNLobjective(model, :Min, objexpr)


    #@NLobjective(m, Max, objexpr)
    println("最適化開始!")
    tic();
    solve(model)
    toc()

    nodes=construct_ancestry(x_t);
    z=getvalue(y)
    println("プログラム完了!")
    println("_____実験結果______")
    println("ASCs")
    println("___________________")
    println("|-----|","True ","|","Est.|")
    for i=1:m
        println("|ASC_",i,"|",b_t[i],"  |",round(z[i],2),"|")
    end
    println("___________________")
    println("Betas")
    println("___________________")
    println("|-----|","True ","|","Est.|")

    for i=m+1:2*m
        println("|β_",i-m,"|",b_t[i],"  |",round(z[i],2),"|")
    end

    println("___________________")
    println("Nest Scales")
    println("___________________")
    println("|---|","True ","|","Est.|")
    for i=1:m-1
        if nodes[i]>0
            println("|μ_",i,"|",mu_t[i],"  |",round(z[(2*m)+i],1),"|")
        end
    end
    println("___________________")
    println("Log-likelihood")
    println("___________________")
    println("|---|","True ","|","Est.|")
    println("|LL""|",round(f_x(z_t),2),"  |",round(likelihood_2(z),2),"|")
    z=getvalue(y)
    l=getobjectivevalue(model)
    return z,l
end

Jacobian_f=z->ForwardDiff.gradient(f,z)
Hessian_f=z->ForwardDiff.hessian(f,z)
Jacobian_f_x=z->ForwardDiff.gradient(f_x,z)



function Jacobian_analytic(z,S,A,X)
    dim1=(1+(m-2));
    dim2=(1+(m-2)+m);
    index_shift=dim1*dim2;
    grad=zeros(length(z))
    grad[1:(end-index_shift)]=Jacobian_f(z[1:(end-index_shift)])
    grad[(end-index_shift+1):end]=gradient_x(z[1:2*m],z[(2*m+1):(end-index_shift)],reshape(z[(end-index_shift+1):end],dim1,dim2),S,A,X)[:]
    return grad
end


function Hessian_analytic(z,S,A,X)
    dim1=(1+(m-2));
    dim2=(1+(m-2)+m);
    index_shift=dim1*dim2;
    l=length(z)
    hess=zeros(l,l);

    y=z[1:(end-index_shift)]
    ###########################ASSIGN BLOCK 1#######################################
    hess[1:(end-index_shift),1:(end-index_shift)]=round.(Hessian_f(y),3);
    ###############################################################################



    ###########################ASSIGN BLOCK 2,3#######################################
    temp_mat=zeros(dim1,dim2,l-index_shift)

    for i=1:1
        for j=2:dim1
            str=string("fd_gradient_rb_",j)
            temp_mat[i,j,:]=feval(str,y)
        end
        for j=dim1+1:dim2
            str=string("fd_gradient_ra_",j-1-(m-2))
            temp_mat[i,j,:]=feval(str,y)
        end
    end

    for j=2:dim1
        for k=2:dim1
            if j==k
                continue
            else
                str=string("fd_gradient_bb_",j,"_",k)
                temp_mat[j,k,:]=feval(str,y)
            end
        end
    end

    for i=2:dim1
        for j=1:m
            str=string("fd_gradient_ba_",i,"_",j)
            temp_mat[i,j+1+(m-2),:]=feval(str,y)
        end
    end

    temp=Array{Float64}(length(y),0)


    for j=1:dim2
        for i=1:dim1
            temp=hcat(temp,temp_mat[i,j,:])
        end
    end

    temp[isnan.(temp)]=0;

    hess[1:(end-index_shift),(end-index_shift)+1:end]=temp
    hess[(end-index_shift)+1:end,1:(end-index_shift)]=temp'

    ########################################################################

    #############################ASSIGN FINAL BLOCK############################
    temp=zeros(dim1,dim2,dim1,dim2)
    b=y[1:p]
    mu=y[p+1:end]

    #ra ra
    for j1=1:m
        for j2=1:m
            temp[1,j1+dim1,1,j2+dim1]=gradient_ra_ra(j1,j2,b,mu,x_t,S_t,A,X)
        end
    end

    #rb ra
    for k=2:dim1
        for j1=1:m
            temp[1,k,1,j1+dim1]=gradient_rb_ra(k,j1,b,mu,x_t,S_t,A,X)
        end
    end

    #rb rb
    for k1=2:dim1
        for k2=2:dim1
            temp[1,k1,1,k2]=gradient_rb_rb(k1,k2,b,mu,x_t,S_t,A,X)
        end
    end

    #rb bb
    for k1=2:dim1
        for l1=2:dim1
            for k2=2:dim1
                temp[k1,l1,1,k2]=gradient_bb_rb(k1,l1,k2,b,mu,x_t,S_t,A,X)
            end
        end
    end

    #rb ba
    for k1=2:dim1
        for j1=1:m
            for k2=2:dim1
                temp[k1,j1+dim1,1,k2]=gradient_ba_rb(k1,j1,k2,b,mu,x_t,S_t,A,X)
            end
        end
    end

    #bb bb
    for k1=2:dim1
        for l1=2:dim1
            if k1==l1
                continue
            end
            for k2=2:dim1
                for l2=2:dim1
                    if k2==l2
                        continue
                    end
                    temp[k1,l1,k2,l2]=gradient_bb_bb(k1,l1,k2,l2,b,mu,x_t,S_t,A,X)
                end
            end
        end
    end

    #gradient_bb_ba(k1,l1,k2,j1,b,mu,x,S,A,X)

    #bb ba
    for k1=2:dim1
        for l1=2:dim1
            if k1==l1
                continue
            end
            for k2=2:dim1
                for j1=1:m
                    temp[k1,l1,k2,j1+dim1]=gradient_bb_ba(k1,l1,k2,j1,b,mu,x_t,S_t,A,X)
                end
            end
        end
    end

    #ba ba
    for k1=2:dim1
        for j1=1:m
            temp[k1,j1+dim1,k1,j1+dim1]=gradient_ba_ba(k1,j1,k1,j1,b,mu,x_t,S_t,A,X)
        end
    end

    hess_temp=Array{Float64}(index_shift,0)

    for j=1:dim2
        for i=1:dim1
            hess_temp=hcat(hess_temp,temp[i,j,:,:][:])
        end
    end

    hess_temp=(hess_temp+hess_temp')-eye(index_shift).*hess_temp
    ##########################################################################
    hess[(end-index_shift)+1:end,(end-index_shift)+1:end]=hess_temp

    hess[isnan.(hess)]=0;



    return hess
end

#Hessian_f=z->ForwardDiff.hessian(f_x,z)
#tic()
#Jacobian_f(z[1])
#toc()
#tic()
#H=Hessian_f(z[1])
#toc()

###############################################################################################
#m=8
#x_t=[
#    0 1 1 1 1 0 0 0 0 0 0 0 0 0 0
#    0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
#    0 0 0 0 0 0 0 0 0 1 1 0 0 0 0
#    0 0 0 0 0 0 0 0 0 0 0 1 1 0 0
#    0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
#    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#    ]
#b_t=[0 1 1 1 1 1 1 1 1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1]'
#mu_t=[1 4 4 4 4 4 4 4 4]

m=8
x_t=[
    0 1 1 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 1 1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 1 1 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
    ]

b_t=[0 1 1 1 1 1 1 1 0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1]'
mu_t=[1 1.5 1.5 2 2 2 2]

#x_t=[
#    0 1 1 0 0 0 0
#    0 0 0 1 1 0 0
#    0 0 0 0 0 1 1
#    ]

#b_t=[0 0.5 0.25 0.3 0.5 -0.1 -0.1 -0.1]'
#mu_t=[1 2 5]


S_t=zeros(m,2*m)
for i=1:m
    S_t[i,i]=1
    S_t[i,i+m]=1
end

data=zeros(n,m)
for i=1:n
    for j=1:m
    d=Normal(0,rand(1)[1]*3)
    data[i,j]=rand(d,1)[1]
    end
end

for i=1:m
    mean_i=mean(data[:,i])
    var_i=var(data[:,i])
    data[:,i]=(data[:,i]-mean_i)./sqrt(var_i)
end


X=[ones(n,m) data];
A=ones(n,m);
R=rand(n,m);

for i=1:n
    for j=1:m
        if R[i,j]>0.9
            A[i,j]=0
        end
    end
end
non_zero_idx=(sum(A,2).>0)[:,1]
A=A[non_zero_idx,:];
X=X[non_zero_idx,:];

n=size(X,1);
dim1=(1+(m-2));
dim2=(1+(m-2)+m);
index_shift=dim1*dim2
x_0=zeros(dim1,dim2);

for i=1:m
    x_0[1,i+1+(m-2)]=1
end



Q=construct_probability_matrix(S_t,X,A,x_t,b_t,mu_t)
C=construct_choice_matrix(Q)
C_true=C

#z_t=hcat(b_t',mu_t,x_t[:]')
#Hessian_analytic(z_t,S_t,A,X)

##############################################################################33
#a,b=solve_NL_subproblem()



#likelihood_x(b_t,mu_t,x_0)

#likelihood_x(b_t,mu_t,x_t)




###############################################################################






function split_vector(x)
    return x[1:2*m],x[(2*m+1):(end-index_shift)],reshape(x[(end-index_shift+1):end],dim1,dim2)
end


function determine_equivalence(x1,x2)
    #println("Testing structure ")
    #print_edges(x1)
    #println("Against structure ")
    #print_edges(x2)
    a1=construct_ancestry(reshape(x1,dim1,dim2))
    a2=construct_ancestry(reshape(x2,dim1,dim2))
    P=collect(1:dim1)

    flag=false

    for j=1+(m-2)+1:1+(m-2)+m
        node1=a1[j]
        node2=a2[j]
        if node1==node2
            continue
        else
            i1=find(y-> (y==node2),a2)
            i2=find(y-> (y==node1),a2)
            a2[i1]=node1
            a2[i2]=node2
            temp=a2[node2]
            a2[node2]=a2[node1]
            a2[node1]=temp

            temp=P[node2]
            P[node2]=P[node1]
            P[node1]=temp
            #push!(nodes_visited,node1)
        end
    end

    if a1==a2
        flag=true
        #println("Equivlalent!!")
    end
    #println("--------------------------------------")
    return flag,P
end


x_true=x_t





###################################################3333
iter_count=0

function solveMasterProblem(T,f_list,x_list,z,H_i,depth,nests,flag)

    println("Solving Master Problem")
    #MasterProblem= Model(solver=KnitroSolver())
    #MasterProblem=Model(solver=CplexSolver(CPXPARAM_OptimalityTarget=3))
    MasterProblem = Model(solver=GurobiSolver())

    @variable(MasterProblem, x[1:dim1,1:dim2], Bin)
    @variable(MasterProblem, μ[1:dim1]>=1)
    @variable(MasterProblem, y[1:m-2],Bin)
    @variable(MasterProblem, β[1:2*m])
    @variable(MasterProblem, η)
    @variable(MasterProblem, δ, Bin)

    #@constraint(MasterProblem,η>=LB+0.01)
    if flag
        @constraint(MasterProblem,η>=0)
    end

    @constraint(MasterProblem,sum(y[i] for i=1:dim1-1)==nests)

    println("Set depth toogled")
    longest_branch = zero(AffExpr)
    for i=1:dim1-1
        longest_branch+=x[i,i+1]
    end
    @constraint(MasterProblem,longest_branch==depth)

    #y=ones(1,m-2)

    #No connections to root
    for i=1:dim1
        @constraint(MasterProblem,x[i,1]==0)
    end

    #No self-arcs
    for i=1:dim1
        @constraint(MasterProblem,x[i,i]==0)
    end

    #Choice nodes belong to one and only one nest
    for i=1:m
        @constraint(MasterProblem,sum(x[j,i+1+(m-2)] for j=1:dim1)==1)
    end


    #Nest nodes can have at most one parent
    for i=1:m-2
        #@constraint(MasterProblem,sum(x[j,i+1] for j=1:dim1)<=1)
        @constraint(MasterProblem,sum(x[j,i+1] for j=1:dim1)==y[i])
        #println(sum(getvalue(x)[j,i+1] for j=1:dim1),getvalue(y)[i])
    end


    #No connections to nest unless chosen
    for i=1:m-2
        for j=1:m-2
            if i==j
                continue
            end
            @constraint(MasterProblem, x[i+1,j+1]+x[j+1,i+1]<=y[i])
            @constraint(MasterProblem, μ[j+1]>=delta+μ[i+1]-mu_upper*(1-x[i+1,j+1]))
        end
        @constraint(MasterProblem, μ[i]<=mu_upper)
    end

    for i=1:m-2
        @constraint(MasterProblem,2*y[i]<=sum(x[i+1,j] for j=1:dim2))
        @constraint(MasterProblem,sum(x[i+1,j] for j=1:dim2)<=m*y[i])
    end
    @constraint(MasterProblem,2<=sum(x[1,j] for j=1:dim2))
    @constraint(MasterProblem,1-dim2*(1-δ)<=sum(x[1,j] for j=1:dim1))
    @constraint(MasterProblem,(m-2)*δ>=sum(y[i] for i=1:m-2))

    for i=1:m
        @constraint(MasterProblem,β[i]<=10)
        @constraint(MasterProblem,β[i]>=-10)
        @constraint(MasterProblem,β[m+i]<=1)
        @constraint(MasterProblem,β[m+i]>=-1)
    end


    @constraint(MasterProblem,β[m]==0)
    @constraint(MasterProblem,μ[1]==mu_norm)

    #Graph is a forest
    @constraint(MasterProblem,sum(sum(x[i,j] for i=1:dim1) for j=1:dim2)==(sum(y[i] for i=1:m-2)+m+1)-1)
    #@constraint(MasterProblem,sum(sum(x[i,j] for i=1:dim1) for j=1:dim2)<=(sum(y[i] for i=1:m-2)+m+1)-1+epsilon)
    #@constraint(MasterProblem,sum(sum(x[i,j] for i=1:dim1) for j=1:dim2)>=(sum(y[i] for i=1:m-2)+m+1)-1-epsilon)
    if length(H_i)>0
        @objective(MasterProblem , Min,η +0.5*(vcat(β,μ,x[:])-z_i)'*(H_i)*(vcat(β,μ,x[:])-z_i))
    else
        @objective(MasterProblem , Min, η)
    end


    #temp=(vcat(β,μ,x[:])-z[1])
    #for j=1:200
        #println(T[1][j]*temp[j])
    #end
    if false
        for i=1:length(T)
            println("i = ",i)
            println(getvalue(η),">=",(vecdot(T[i],(vcat(getvalue(β),getvalue(μ),getvalue(x)[:])-z[i])),"+",f_x(z[i])))
            println(getvalue(η),"<=",f_x(z[i])-1)
            println("---------------------------------------------------------")
        end
    end


    for i=1:length(T)
        @constraint(MasterProblem, η>=(vecdot(T[i],(vcat(β,μ,x[:])-z[i]))+f_list[i]))
        @constraint(MasterProblem, η<=f_list[i]-0.01)
        x_temp=reshape(x_list[i],dim1,dim2)
        @constraint(MasterProblem,sum(sum((2*x_temp[j,k]-1)*x[j,k] for j=1:dim1) for k=1:dim2) <= sum(sum((x_temp[j,k]) for j=1:dim1) for k=1:dim2)-1)
    end

    b=0
    if b>0
        for j=1:m
            @constraint(MasterProblem,β[j]==z[b][j])
            @constraint(MasterProblem,β[m+j]==z[b][m+j])
        end
        for j=1:1+m-2
            @constraint(MasterProblem, μ[j]==z[b][j+2*m])
        end
    end

    global iter_count=0

    function cycle_detection_2(model)
        println("----\nInside cycle detection 2.0 callback")
        A=convert(Array{Int64},vcat(round.(getvalue(x),2),zeros(10,19)))
        Tree=DiGraph(A)
        if iter_count>500
            println("Exceeded count :(")
            return
        end

        components=weakly_connected_components(Tree)

        if length(components) != m-2 - sum(getvalue(y))+1
            print("Found disconnected components")
            global iter_count+=1

            for component in components
                if length(component)>1
                    arcs_in_component = zero(AffExpr)
                    nodes_in_component= zero(AffExpr)

                    for node in component
                        if node <= 1+(m-2)
                            if node==1
                                nodes_in_component+=1
                            else
                                nodes_in_component+=y[node-1]
                            end
                        else
                            nodes_in_component+=1
                        end
                        node_neighborhood=neighborhood(Tree,node,1)
                        if length(node_neighborhood)>1
                            for neighbor in node_neighborhood
                                arcs_in_component+=x[node,neighbor]
                            end
                        end
                    end

                    @lazyconstraint(model, arcs_in_component <= nodes_in_component-1)
                end

            end
        else
            print("Graph is connected")
            return
        end

    end


    function cycle_detection(model)
        println("----\nInside cycle detection callback")
        #A=convert(Array{Int64},vcat(getvalue(x),zeros(10,19)))
        #A=convert(Array{Int64},vcat(floor.(getvalue(x)),zeros(10,19)))
        c_x=round.(getvalue(x),1)
        A=convert(Array{Int64},vcat(c_x,zeros(m,dim2)))
        #A=vcat(getvalue(x),zeros(10,19))
        Tree=DiGraph(A.|A')
        #Tree=DiGraph(A)

        if iter_count>500
            println("Exceeded count :(")
            return
        end

        if is_cyclic(Tree)
            global iter_count+=1

            cycles= simplecycles(Tree)
            for cycle in cycles
                if length(cycle)<3
                    continue
                else
                    arcs_in_cycle = zero(AffExpr)
                    nodes_in_cycle= zero(AffExpr)
                    cycle_length=length(cycle)
                    println("A cycle of length ",cycle_length," has been detected")

                    for i=1:cycle_length
                        #println("cycle ",i, "is ", cycle[i])
                        nodes_in_cycle+=y[cycle[i]-1]
                    end

                    arcs_in_cycle+=c_x[cycle[1],cycle[2]]*x[cycle[1],cycle[2]]+c_x[cycle[2],cycle[1]]*x[cycle[2],cycle[1]]
                    for i=2:cycle_length-1
                        arcs_in_cycle+=c_x[cycle[i],cycle[i+1]]*x[cycle[i],cycle[i+1]]+c_x[cycle[i+1],cycle[i]]*x[cycle[i+1],cycle[i]]
                    end
                    arcs_in_cycle+=c_x[cycle[cycle_length],cycle[1]]*x[cycle[cycle_length],cycle[1]]+c_x[cycle[1],cycle[cycle_length]]*x[cycle[1],cycle[cycle_length]]
                    println("Adding cycle elimination cut")
                    @lazyconstraint(model, arcs_in_cycle <= nodes_in_cycle-1)

                end

            end
        else
            println("No cycles detected!")
            return
        end

    end

    function equivalence_cut(model)
        println("----\nInside Equivalence cut generator")

        c_x=round.(getvalue(x),1)

        flag=false
        for i=1:length(T)
            flag, P=determine_equivalence(c_x[:],x_list[i])
            if flag==true
                println("Equivalent structure found!")
                println("Current")
                print_edges(c_x[:])
                println("------------")

                println("Previously estimated")
                print_edges(x_list[i])
                #beta_temp, mu_temp, x_temp=split_vector(z[i])
                #dbeta_temp, dmu_temp, dx_temp=split_vector(T[i])

                #mu_temp=mu_temp[P]
                #x_temp=x_temp[P,vcat(P,collect(dim1+1:dim2))]

                #dmu_temp=dmu_temp[P]
                #dx_temp=dx_temp[P,vcat(P,collect(dim1+1:dim2))]

                #if c_x!=x_temp
                    #println("MONDAI !!")
                #end

                #println("Adding Equivalence cuts")
                #z_temp=vcat(beta_temp,mu_temp,c_x[:])
                #gradient_temp=Jacobian_f_x(z_temp)
                #gradient_temp=Jacobian_analytic(z_temp,S_t,A,X)
                #f=f_x(z_temp)
                #gradient_temp[isnan.(gradient_temp)]=0
                @lazyconstraint(model,sum(sum((2*c_x[j,k]-1)*x[j,k] for j=1:dim1) for k=1:dim2) <= sum(sum((c_x[j,k]) for j=1:dim1) for k=1:dim2)-1)
                #@lazyconstraint(model,η>=(vecdot(gradient_temp,(vcat(β,μ,x[:])-z_temp))+f))
                #@lazyconstraint(model,η<=f-0.01)
            end
        end
        println("No more cuts!")
        return
    end

    function depth_cut(model)
        println("Inside pruning callback")
        current_x=round.(getvalue(x),1)
        #println(current_x)
        #print_edges(current_x)
        nodes=construct_ancestry(current_x)
        for j=dim1+1:dim2
            #println("looking at alternative ",j)
            path=path_from_root(j,nodes)
            #println("path computed")
            l=length(path)
            depth_path=l-2
            if depth_path>depth
                println("Maximum depth Exceeded")

                #@lazyconstraint(model,sum(sum((2*c_x[j,k]-1)*x[j,k] for j=1:dim1) for k=1:dim2) <= sum(sum((c_x[j,k]) for j=1:dim1) for k=1:dim2)-1)

                branch = zero(AffExpr)
                for i=1:l-2
                    branch+=x[path[i],path[i+1]]
                end
                println("Adding prunning constraint")
                @lazyconstraint(model,branch<=depth)
            end
        end
        println("Tree pruned!")
        return
    end
    addlazycallback(MasterProblem,cycle_detection)
    addlazycallback(MasterProblem,equivalence_cut)
    addlazycallback(MasterProblem,depth_cut)
    status=solve(MasterProblem)
    return status, getvalue(x), getvalue(η)
end

#z_0=hcat(b_t',mu_t,x_0[:]')
z_t=hcat(b_t',mu_t,x_true[:]')
#status,start_x, temp_eta=solveMasterProblem([],[],[],1,4,true)

function optimal_logit(depth,nests,max_iter,flag)
    count=0
    #get initial solution
    status,start_x, temp_eta=solveMasterProblem([],[],[],[],[],depth,nests,true)
    if status!=:Optimal
        return 0,start_x,count
    end

    global x_t=start_x

    f_list=[]
    η_list=[]
    x_list=[]
    T=[];
    z=[];
    push!(x_list,x_t[:]')

    for i=2:max_iter
        count+=1
        temp, f_i=solve_NL_subproblem()
        #push!(WARM_START,temp)
        beta_i=temp[1:2*m]
        mu_i=temp[(2*m+1):end]
        #z_i=hcat(b_t',mu_t,x_0[:]')'
        #f_i=10000
        z_i=hcat(beta_i',mu_i',x_t[:]')'
        push!(f_list,f_i)
        push!(z,z_i)

        println("Calculating gradient")
        push!(T,Jacobian_analytic(z_i,S_t,A,X))
        println("Done")
        #push!(T,Jacobian_f_x(z_i))

        #T[i][1:2*m+1+(m-2)]=T[i][1:2*m+1+(m-2)]/100

        #println("Calculating Hessian")
        #H_i=Hessian_analytic(z_i,S_t,A,X)
        status, x_i,η_i=solveMasterProblem(T,f_list,x_list,z,[],depth,nests,false)
        if status!=:Optimal
            print("Infeasibility")
            break
        end
        x_i=round.(x_i,1)
        println("Current η is: ",η_i)


        push!(η_list,η_i)


        push!(x_list,x_i[:]')
        println("Printing tree @ iteration ",i)
        print_edges(x_list[i])
        println("-------------------------------")
        global x_t=x_i
    end
    a,b=findmin(f_list)
    if flag
        return f_list,η_list,count
        #z_t=hcat(b_t',mu_t,x_true[:]')
        plotly()
        #x_t=x_true
        #temp, f_i=solve_NL_subproblem()
        plot(f_list,xlabel="Iteration", ylabel="Objective Value",lab="NLP SubProblem")
        plot!(η_list, lab="MILP Master Problem")
        #plot!(f_i*ones(iter_num,1),lab="Optimal Value")
    end
    return f_list[b],z[b],count
end



nest_list=0:1:m-2
depth_list=0:1:m-2
max_iter=15
p=length(nest_list)
q=length(depth_list)
total_iterations=0
LL = zeros(p, q)
C_X=[]
for i in 1:p
    for j in 1:q
        if depth_list[j]>nest_list[i]
            continue
        end
        println("-----------------------------------------------")
        println("Nests = ",nest_list[i], " Depth = ",depth_list[j])
        println("-----------------------------------------------")
        optimal_value,optimal_solution, itr_count=optimal_logit(depth_list[j],nest_list[i],max_iter,false)
        LL[i, j] = optimal_value
        push!(C_X,optimal_solution)
        total_iterations+=itr_count
    end
end
temp_x=[]
count=0
for i in 1:p
    for j in 1:q
        if depth_list[j]>nest_list[i]
            continue
        end
        count+=1
        if LL[i,j]<712 && LL[i,j]>0
            println("-----------------------------------------------")
            println("Nests = ",nest_list[i], " Depth = ",depth_list[j])
            println("-----------------------------------------------")
            println(LL[i,j])
            temp_x=C_X[count]
        end

    end
end

for i in 1:p
    for j in 1:q
        if LL[i,j]==0
            LL[i,j]=NaN
        end
    end
end

x=[]
y=[]
z=[]

for i in 1:p
    for j in 1:q
        if depth_list[j]>nest_list[i]
            continue
        end
        push!(x,nest_list[i])
        push!(y,depth_list[j]+1)
        push!(z,LL[i,j])
    end
end




plotly()
xgrid, ygrid = meshgrid(nest_list, depth_list)
Plots.plot([x y z],xlabel="Nests", ylabel="Depth",zlabel="ℒ")
Plots.scatter(x,y,z,xlabel="Nests", ylabel="Depth",zlabel="ℒ")

#Plots.surface!(xgrid, ygrid, z2', lab="Optimal value")

heatmap(nest_list, depth_list+1, LL', aspect_ratio=1,xlabel="Number of Nests", ylabel="Tree Depth",zlabel="ℒ", title="-ℒ (training)")


X_true=X

data=zeros(n,m)
for i=1:n
    for j=1:m
    d=Normal(0,rand(1)[1]*3)
    data[i,j]=rand(d,1)[1]
    end
end

for i=1:m
    mean_i=mean(data[:,i])
    var_i=var(data[:,i])
    data[:,i]=(data[:,i]-mean_i)./sqrt(var_i)
end

X=[ones(n,m) data];
A=ones(n,m);
R=rand(n,m);

for i=1:n
    for j=1:m
        if R[i,j]>0.9
            A[i,j]=0
        end
    end
end
Q_test=construct_probability_matrix(S_t,X,A,x_true,b_t,mu_t)
C_test=construct_choice_matrix(Q_test)

C=C_test

x_test=[]
y_test=[]
LL_test=zeros(p, q)

count=0
for i in 1:p
    for j in 1:q
        if depth_list[j]>nest_list[i]
            continue
        end
        count+=1
        if isnan(C_X[count][1])
            continue
        end


        LL_test[i,j]=f_x(C_X[count])
    end
end


for i in 1:p
    for j in 1:q
        if LL_test[i,j]==0
            LL_test[i,j]=NaN
        end
    end
end

plotly()
upscale = 1 #8x upscaling in resolution
fntsm = Plots.font("sans-serif", 8*upscale)
fntlg = Plots.font("sans-serif", 8*upscale)
#Plots.scalefontsizes(8)
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
default(size=(800*upscale,600*upscale)) #Plot canvas size
default(dpi=300) #Only for PyPlot - presently broken

heatmap(nest_list, depth_list+1, LL_test', aspect_ratio=1,xlabel="Number of Nests", ylabel="Tree Depth",zlabel="ℒ", title="-ℒ (test)")

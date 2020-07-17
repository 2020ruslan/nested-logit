const delta=0.0;
const mu_norm=1;
const mu_upper=20;
const epsilon=1e-20;
env = Gurobi.Env()

######################Graph related functions#######################

function construct_ancestry(x)
    dim1,dim2=size(x)
    nodes=zeros(dim2,1);
    for i=2:dim2
        parent=findall(y-> (y==1),x[:,i]);
        if length(parent)>0
            nodes[i]=findall(y-> (y==1),x[:,i])[1];
        else
            nodes[i] =0
        end
    end
    nodes=convert(Array{Int64},nodes)
    return nodes
end


function get_children(nest,x)
    dim1,dim2=size(x)
    return [child for child ∈ 1:dim2 if x[nest,child]>0]
end

function construct_graph_from_ancestry(ancestry)
    dim2=length(ancestry)
    dim1=Int64(0.5*(dim2-1))
    graph=zeros(dim1,dim2)

    for node ∈ 1:dim2
        parent=ancestry[node]
        if parent>0
            graph[parent,node]=1
        end
    end
    return graph
end

function construct_depth_map(x)
    dim1,~=size(x)
    depth_map=zeros(dim1,1);
    nodes=construct_ancestry(x)
    depth_map[1]=0
    for i=2:dim1
        if nodes[i]==0
            depth=-1
        else
            depth=0;
            parent=nodes[i]
            while parent>0
                depth+=1
                parent=nodes[parent]
            end
        end
        depth_map[i]=depth
    end
    return depth_map
end


function depth_sorted_nest_list(x)
    nodes=construct_ancestry(x)
    nest_list=get_nests(1,x)
    depth_map=construct_depth_map(x)
    depth_list=[depth_map[nest] for nest in nest_list]
    sort_idx=sortperm(depth_list,rev=true)
    sorted_nest_list=nest_list[sort_idx]
    return sorted_nest_list,nodes
end


function nest_children_bool(x)
    dim1,dim2=size(x)
    nest_children_map=zeros(dim1,1);
    for i=1:dim1
        nest_children=findall(y-> (y==1),x[i,1:dim1])
        if nest_children==[]
            nest_children_map[i]=0
        else
            nest_children_map[i]=1
        end
    end
    return nest_children_map.>0
end

function print_edges(x)
    dim1,dim2=size(x)
    x_val=(round.(x,digits=1).==1)*1
    #A=convert(Array{Int64},vcat(x_val,zeros(dim1+1,dim2)))
    #A=convert(Array{Int64},vcat(round.(reshape(x,dim1,dim2),digits=1),zeros(dim1+1,dim2)))
    #Tree=DiGraph(A)
    nodes=construct_ancestry(x_val)
    for i∈ 2:dim1
        if nodes[i]>0
            if nodes[i]==1
                println("Root to nest ",i-1)
            else
                println("Nest ",nodes[i]-1," to nest ",i-1)
            end
        end
    end
    for i∈ dim1+1:dim2
        if nodes[i]==1
            println("Root to alternative ",i-dim1)
        else
            println("Nest ",nodes[i]-1," to alternative ",i-dim1)
        end
    end
    #for e in edges(Tree)
     #  println(e)
    #end
end

function draw_tree(x,mu)
    dim1,dim2=size(x)
    m=dim1+1
    nodes=construct_ancestry(x)
    Tree=DiGraph(convert(Array{Int64},vcat(x,zeros(m,dim2))))
    labels=["r"]
    for nest ∈ 2:dim1
        if nodes[nest]>0
            push!(labels,string("μ= ",mu[nest]))
        else
            push!(labels,"")
        end
    end
    for alternative ∈ 1:m
        push!(labels,string(alternative))
    end
    TikzGraphs.plot(Tree,labels, node_style="draw,rounded corners, fill=white!10",edge_style="white",options="scale=2, font=\\huge")

end


function path_from_root(j,nodes)
    path_j=[]
    node=j
    push!(path_j,node)
    while node>0
        node=nodes[node]
        push!(path_j,node)
        if length(path_j)>=num_alternatives
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


function get_nests(nest,x)
    dim1,~=size(x)
    nests=[]
    if nest==1
        append!(nests,1)
    end
    progeny=findall(y-> (y==1),x[nest,1:dim1])
    append!(nests,progeny)
    for child in progeny
        append!(nests,get_nests(child,x))
    end
    return nests
end
function is_valid_tree(x)
    dim1,dim2=size(x)
    ancestry=construct_ancestry(x)

    if sum(x)!=sum(sum(x[i,:])>0 for i=1:dim1)+dim1
        return false
    end

    for child ∈ 2:dim2
        if sum(x[parent,child] for parent ∈ 1:dim1)>1
            return false
        end
    end

    for child ∈ dim1+1:dim2
        if ancestry[child]==0
            return false
        end

    end
    return true
end

function determine_equivalence(x1,x2)
    dim1,dim2=size(x1)

    a1=construct_ancestry(x1)
    a2=construct_ancestry(x2)
    P=collect(1:dim1)

    #println("a1 ", a1)
    #println("x1", x1)
    #println("a2 ", a2)
    #println("x2", x2)

    if !(is_valid_tree(x1) && is_valid_tree(x2))
        #println("Invalid tree passed")
        return true
    end

    flag=false

    for j=dim1+1:dim2
        parent_of_j_in_1=a1[j]
        parent_of_j_in_2=a2[j]

        if parent_of_j_in_1==parent_of_j_in_2
            continue
        else
            siblings_of_j_in_2= a2.==parent_of_j_in_2
            children_of_parent_of_j_in_1_in_2= a2.==parent_of_j_in_1

            a2[siblings_of_j_in_2].=parent_of_j_in_1
            a2[children_of_parent_of_j_in_1_in_2].=parent_of_j_in_2

            temp=a2[parent_of_j_in_2]
            a2[parent_of_j_in_2]=a2[parent_of_j_in_1]
            a2[parent_of_j_in_1]=temp

            temp=P[parent_of_j_in_2]
            P[parent_of_j_in_2]=P[parent_of_j_in_1]
            P[parent_of_j_in_1]=temp
            #push!(nodes_visited,node1)
        end
    end

    if a1==a2
        flag=true
    end
    return flag
end


function split_vector(x)
    dim1,dim2=size(x)
    index_shift=dim1*dim2
    return x[1:2*m],x[(2*m+1):(end-index_shift)],reshape(x[(end-index_shift+1):end],dim1,dim2)
end
##############################################################################


########################################Data preperation functions###############

function split_at(at,n)
    idx = shuffle(1:n)
    train_idx = view(idx, 1:floor(Int, at*n))
    test_idx = view(idx, (floor(Int, at*n)+1):n)
    return train_idx, test_idx
end

function generate_montecarlo_data(n,S,x,b,mu,γ,train_percent)
    m,p=size(S)
    @assert p>=m "Misspecified Model"

    data=zeros(n,p-m)

    @simd for j=1:p-m
        d=Normal(0,1)
        data[:,j]=rand(d,n)

    end

    X=[ones(n,m) data];
    A=ones(n,m);
    R=rand(n,m);

    for i=1:n
        for j=1:m
            if R[i,j]>γ
                A[i,j]=0
            end
        end
    end
    non_zero_idx=(sum(A,dims=2).>0)[:,1]
    A=A[non_zero_idx,:];
    X=X[non_zero_idx,:];


    Q=construct_probability_matrix(S,X,A,x,b,mu)
    C,Y=construct_choice_matrix(Q)

    #histogram(Y,title="Choices")
    #histogram(A,title= "Availabilities")

    train_idx, test_idx=split_at(train_percent,sum(non_zero_idx))

    A_train=A[train_idx,:]
    C_train=C[train_idx,:]
    X_train=X[train_idx,:]
    Y_train=Y[train_idx,:]

    A_test=A[test_idx,:]
    C_test=C[test_idx,:]
    X_test=X[test_idx,:]
    Y_test=Y[test_idx,:]

    return A_train,C_train,X_train,Y_train,A_test,C_test,X_test,Y_test

end


#######################################Likelihood related functions###########
#function V(alternative,S,X,b,p)
#    if alternative==1
#        return b[1] + b[7]*X[:DIST]/3.1
#    elseif alternative==2
#        return b[2] + b[8]*X[:DIST]/9.6 + b[13]*X[:CBD] +b[18]*X[:TRANS]
#    elseif alternative==3
#        return b[3] + b[9]*X[:CBD] + b[14]*X[:TRANS]  + b[19]*X[:CAR_COST] + b[23]*X[:CAR_TT]
#    elseif alternative==4
#        return b[4] + b[10]*X[:CBD] + b[15]*X[:TRANS] + b[20]*X[:CAR_COST2] + b[24]*X[:CAR_TT]
#    elseif alternative==5
#        return b[5] + b[11]*X[:CBD] +  b[16]*X[:TRANS] + b[21]*X[:CAR_COST3] + b[25]*X[:CAR_TT]
#    elseif alternative==6
#        return b[6] + b[12]*X[:CBD] +  b[17]*X[:TRANS]+ b[22]*X[:PT_COST] + b[26]*X[:PT_TT]
#    else
#        println("Error: Invalid index for alternative")
#        return 0
#    end

#end

function V(alternative,S,X,b,p)
      return sum(b[j]*S[alternative,j]*X[j] for j=1:p)
end



function Inclusive_value(IV_lookup_table,ind_idx,i,S,X,A,x,b,mu)
    dim1,dim2=size(x)
    IV=epsilon
    for j=1+dim1:dim2
        if A[j-dim1]*x[i,j]>0
            IV+= exp(mu[i]*IV_lookup_table[ind_idx,j])
        end
    end
    for j=2:dim1
        if x[i,j]>0
            IV+= exp(mu[i]*IV_lookup_table[ind_idx,j])
        end
    end
    return log(IV)/mu[i]
end

function IV_AD(nest,ind_idx,S,X,A,x,b,mu,p)
    dim1,dim2=size(x)
    IV=epsilon
    for j=1:m
        if A[ind_idx,j]*x[nest,j]>0
            IV+=exp(mu[nest]*V(j,S,X[ind_idx,:],b,p))
        end
    end
    for j=2:dim1
        if x[nest,j]>0
            IV+=x[nest,j]*exp(mu[nest]*IV_AD(j,ind_idx,S,X,A,x,b,mu,p))
        end
    end
    return log(IV)/mu[nest]
end


function Probability_AD(nodes,individual,alternative,S,X,A,Y,x,b,mu,p,dim1)
    parent=nodes[alternative]
    temp=mu[nodes[alternative]]*V(alternative-dim1,S,X[individual,:],b,p)
    anscestor=nodes[parent]
    while anscestor>0
        temp=temp+(mu[anscestor]-mu[parent])*IV_AD(parent,individual,S,X,A,x,b,mu,p)
        parent=anscestor
        anscestor=nodes[parent]
    end
    temp=temp-mu[parent]*IV_AD(parent,individual,S,X,A,x,b,mu,p)
    return temp
end


function likelihood_AD(S,X,A,Y,x,b,mu)
    n,~=size(A)
    dim1,~=size(x)
    nodes=construct_ancestry(x)
    m,p=size(S)
    ll=0
    for i=1:n
        ll+=Probability_AD(nodes,i,Y[i]+dim1,S,X,A,Y,x,b,mu,p,dim1)
    end
    return ll
end


function full_gradient_AD(b,mu,x,S,A,X,Y)
    dim1,dim2=size(x)
    m,p=size(S)
    index_shift=dim1*dim2;
    function AD_likelihood(y)
        return -1*likelihood_AD(S,X,A,Y,reshape(y[(end-index_shift+1):end],dim1,dim2),y[1:p],y[(p+1):end-index_shift])
    end
    Jacobian_likelihood=z->ForwardDiff.gradient(AD_likelihood,z)

    return reshape(Jacobian_likelihood(hcat(b',mu',x[:]'))[(end-index_shift+1):end],dim1,dim2)
end


function construct_IV_lookup_table(sorted_nest_list,S,X,A,x,b,mu)
    n,m=size(A)
    m,p=size(S)
    dim1=1+(m-2)
    dim2=dim1+m

    IV_lookup_table=Matrix{Any}(undef,n,dim2)
    #nest_availiability=Array{Bool}(undef,n,dim1)
    for i=1:n
        for alternative=1+dim1:dim2
            if A[i,alternative-dim1]>0
                IV_lookup_table[i,alternative]=V(alternative-dim1,S,X[i,:],b,p)
            end
        end
        for nest in sorted_nest_list
            IV_lookup_table[i,nest]=Inclusive_value(IV_lookup_table,i,nest,S,X[i,:],A[i,:],x,b,mu)
            #nest_availiability[i,nest]=any(x[nest,j+dim1]*A[i,j]>0 for j=1:m) || any(x[nest,j]*nest_availiability[i,j]>0 for j=1:dim1)
            #println(nest_availiability[i,nest])
            #if nest_availiability[i,nest]
                #IV_lookup_table[i,nest]=Inclusive_value(IV_lookup_table,i,nest,S,X[i,:],A[i,:],x,b,mu)
            #else
                #IV_lookup_table[i,nest]=0
            #end
        end
    end
    return IV_lookup_table
end


function Probability(IV_lookup_table,nodes,i,j,mu)
    parent=nodes[j]
    temp=mu[nodes[j]]*IV_lookup_table[i,j]
    anscestor=nodes[parent]
    while anscestor>0
        temp=temp+(mu[anscestor]-mu[parent])*IV_lookup_table[i,parent]
        parent=anscestor
        anscestor=nodes[parent]
    end
    temp=temp-mu[parent]*IV_lookup_table[i,parent]
    return temp
end


function construct_probability_matrix(S,X,A,x,b,mu)
    n,m=size(A)
    dim1,dim2=size(x)
    nodes=construct_ancestry(x);
    Q=zeros(n,m);
    sorted_nest_list,nodes=depth_sorted_nest_list(x)
    IV_lookup=construct_IV_lookup_table(sorted_nest_list,S,X,A,x,b,mu)
    for i=1:n
        for j=1+dim1:dim2
            if A[i,j-dim1]>0
                Q[i,j-dim1]=exp(Probability(IV_lookup,nodes,i,j,mu))
            end
        end
    end
    return Q./sum(Q,dims=2)
end


function construct_choice_matrix(Q)
    n,m=size(Q)
    Y=m.-sum(rand(n,1).<cumsum(Q,dims=2),dims=2).+1;
    C=zeros(n,m);
    for i=1:n
        for j=1:m
            if j==Y[i,1]
                C[i,j]=1;
                break
            end
        end
    end
    return C,Y
end

function likelihood(S,X,A,Y,x,sorted_nest_list,nodes,b,mu)
    n,~=size(A)
    dim1,~=size(x)

    IV_lookup=construct_IV_lookup_table(sorted_nest_list,S,X,A,x,b,mu)

    ll=Array{Any}(undef,n)
    Threads.@threads for i=1:n
        ll[i]=Probability(IV_lookup,nodes,i,Y[i]+dim1,mu)
    end
    return sum(ll)
end



function neg_likelihood(S,X,A,Y,x,sorted_nest_list,nodes,b,mu)
    return -1*likelihood(S,X,A,Y,x,sorted_nest_list,nodes,b,mu)
end


###################################Strucutral gradient related functions##########
#root r to alternative a
function gradient_root_to_alternative(a,mu,A,Y,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    grad=0
    @simd for i=1:n
        if A[i,a-dim1]>0
            util=IV_lookup[i,a]
            utilr=IV_lookup[i,1]
            grad+=-exp(mu[1]*(util-utilr))
            if Y[i]+dim1==a
                grad+=mu[1]*(util-utilr)
            end
        end
    end
    return grad
end

function gradient_root_to_alt_root_to_alt(a1,a2,mu,A,Y,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    grad=0
    @simd for i=1:n
        if A[i,a1-dim1]&&A[i,a2-dim1]>0
            utila1=IV_lookup[i,a1]
            utila2=IV_lookup[i,a2]
            utilr=IV_lookup[i,1]
            grad+=+exp(mu[1]*(utila1-utilr))*exp(mu[1]*(utila2-utilr))
            if Y[i]+dim1==a1
                grad+=-exp(mu[1]*(utila2-utilr))
            end
        end
    end
    return grad
end

function gradient_root_to_alt_root_to_nest(a,b,mu,A,Y,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    grad=0
    @simd for i=1:n
        if A[i,a-dim1]>0
            utila=IV_lookup[i,a]
            utilb=IV_lookup[i,b]
            utilr=IV_lookup[i,1]
            grad+=+exp(mu[1]*(utila-utilr))*exp(mu[1]*(utilb-utilr))
            if Y[i]+dim1==a
                grad+=-exp(mu[1]*(utilb-utilr))
            end
        end
    end
    return grad
end


function gradient_root_to_alt_nest_to_alt(a1,b,a2,mu,x,A,Y,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    grad=0
    @simd for i=1:n
        if A[i,a1-dim1]&&A[i,a2-dim1]>0
            utila1=IV_lookup[i,a1]
            derivative=δΓ(1,b,a2,i,mu,x,IV_lookup)
            utilr=IV_lookup[i,1]
            grad+=+mu[1]*exp(mu[1]*(utila1-utilr))*derivative
            if Y[i]+dim1==a
                grad+=-mu[1]*derivative
            end
        end
    end
    return grad
end


function gradient_root_to_alt_nest_to_nest(a,b1,b2,mu,x,A,Y,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    grad=0
    @simd for i=1:n
        if A[i,a-dim1]>0
            utila=IV_lookup[i,a]
            derivative=δΓ(1,b1,b2,i,mu,x,IV_lookup)
            utilr=IV_lookup[i,1]
            grad+=+mu[1]*exp(mu[1]*(utila-utilr))*derivative
            if Y[i]+dim1==a
                grad+=-mu[1]*derivative
            end
        end
    end
    return grad
end

#root r to nest
function gradient_root_to_nest(nest,mu,Y,nodes,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    m=dim1+1

    if nodes[nest]==0
        return 0
    end

    adjusted_ancestry=copy(nodes)
    adjusted_ancestry[nest]=1

    nest_in_alt_ancestry=Array{Bool}(undef,m)
    for a=1:m
        nest_in_alt_ancestry[a],~=path_has_node(path_from_root(a+dim1,nodes),nest)
    end

    grad=0
    @simd for i=1:n
        utilk=IV_lookup[i,nest]
        utilr=IV_lookup[i,1]
        norm_const=utilk-utilr
        grad+=-1*exp(mu[1]*norm_const)
        a=Y[i]

        if nest_in_alt_ancestry[a]
            grad+=Probability(IV_lookup,adjusted_ancestry,i,a+dim1,mu)
        end
    end
    return grad
end

#root r to nest
function gradient_root_to_nest_root_to_nest(nest1,nest2,mu,x,Y,nodes,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    m=dim1+1

    if nodes[nest]==0
        return 0
    end

    adjusted_ancestry=copy(nodes)
    adjusted_ancestry[nest1]=1

    nest_in_alt_ancestry=Array{Bool}(undef,m)
    for a=1:m
        nest_in_alt_ancestry[a],~=path_has_node(path_from_root(a+dim1,nodes),nest1)
    end

    grad=0
    @simd for i=1:n
        utilk=IV_lookup[i,nest1]
        utilr=IV_lookup[i,1]
        derivative=δΓ(1,nest1,nest2,i,mu,x,IV_lookup)
        norm_const=utilk-utilr
        grad+=mu[1]*exp(mu[1]*norm_const)*derivative
        a=Y[i]

        if nest_in_alt_ancestry[a]
            grad+=-1*mu[1]*derivative
        end
    end
    return grad
end


#δΓa wrt x_bc
# a: nest
# b: nest
# c: node
function δΓ(a,b,c,person,μ,x,IV_lookup)
    dim1=length(μ)
    if a==b
        util_diff=IV_lookup[person,c]-IV_lookup[person,a]
        #println("P1")
        return (1/μ[a])*exp(μ[a]*util_diff)
    else
        grad=0
        for k=2:dim1
            #println("P2")
            if x[a,k]>0
                #println("P2")
                util_diff=IV_lookup[person,k]-IV_lookup[person,a]
                grad+=x[a,k]*μ[a]*exp(μ[a]*util_diff)*δΓ(k,b,c,person,μ,x,IV_lookup)
            end
        end
        return (1/μ[a])*grad
    end
end

function accumulate_δΓ(nodes,derivative_nest,path_alternative,derivative_alternative,person,mu,x,IV_lookup)
    temp=0
    parent=nodes[path_alternative]
    anscestor=nodes[parent]
    while anscestor>0
        temp=temp+(mu[anscestor]-mu[parent])*δΓ(parent,derivative_nest,derivative_alternative,person,mu,x,IV_lookup)
        parent=anscestor
        anscestor=nodes[parent]
    end
    temp=temp-mu[1]*δΓ(1,derivative_nest,derivative_alternative,person,mu,x,IV_lookup)
    return temp
end

function gradient_nest_to_alternative(nest,alternative,mu,A,Y,nodes,x,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)

    if nodes[nest]==0
        return 0
    end

    adjusted_ancestry=copy(nodes)
    adjusted_ancestry[alternative]=nest

    #flag,~=path_has_node(path_from_root(alternative,nodes),nest)

    grad=0
    @simd for i=1:n
        if A[i,alternative-dim1]>0
            grad+=accumulate_δΓ(nodes,nest,alternative,alternative,i,mu,x,IV_lookup)
            a=Y[i]
            if alternative==a+dim1
                grad+=Probability(IV_lookup,adjusted_ancestry,i,alternative,mu)
            end
        end
    end

    return grad
end

function gradient_root_to_nest_to_alternative(rnest,nest,alternative,mu,A,Y,nodes,x,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)

    if nodes[nest]==0
        return 0
    end

    adjusted_ancestry=copy(nodes)
    adjusted_ancestry[alternative]=nest

    #flag,~=path_has_node(path_from_root(alternative,nodes),nest)
    grad=0
    @simd for i=1:n
        if A[i,alternative-dim1]>0
            utilr=IV_lookup[i,1]
            utilrnest=IV_lookup[i,rnest]
            deriv_diff=δΓ(rnest,nest,alternative,i,mu,x,IV_lookup)-δΓ(1,nest,alternative,i,mu,x,IV_lookup)
            grad+=-mu[1]*exp(mu[1]*(utilrnest-utilr))*deriv_diff

            a=Y[i]
            if alternative==a+dim1
                grad+=-exp(mu[1]*(utilrnest-utilr))
                if x[1,rnest]>0
                    grad+=Probability(IV_lookup,adjusted_ancestry,i,alternative,mu)
                end
            end
        end
    end

    return grad
end

function gradient_nest_to_nest(nest_origin,nest_destination,mu,Y,nodes,x,IV_lookup)
    dim1=length(mu)
    n,~=size(Y)
    m=dim1+1

    #println("flag:", nodes)


    #println(nodes[nest_origin]==0 || nodes[nest_destination]==0)

    if nodes[nest_origin]==0 || nodes[nest_destination]==0
        return 0
    end

    adjusted_ancestry=copy(nodes)
    adjusted_ancestry[nest_destination]=nest_origin

    #dest_nest_in_alt_ancestry=Array{Bool}(undef,m)
    #for a=1:m
        #dest_nest_in_alt_ancestry[a],~=path_has_node(path_from_root(a+dim1,nodes),nest_destination)
    #end

    #path_has_node(path_from_root(nest_destination,nodes),nest_origin)

    flag,~=path_has_node(path_from_root(nest_origin,nodes),nest_destination)


    nest_in_alt_ancestry=Array{Bool}(undef,m,2)
    nest_position_in_alt_ancestry=zeros(m,2)

    for a=1:m
        nest_in_alt_ancestry[a,1],nest_position_in_alt_ancestry[a,1]=path_has_node(path_from_root(a+dim1,nodes),nest_origin)
        nest_in_alt_ancestry[a,2],nest_position_in_alt_ancestry[a,2]=path_has_node(path_from_root(a+dim1,nodes),nest_destination)
    end

    grad=0
    @simd for i=1:n
        a=Y[i]
        #println(nodes,"-",nest_origin,"-",a+dim1,"-",nest_destination,"-",i)
        grad+=accumulate_δΓ(nodes,nest_origin,a+dim1,nest_destination,i,mu,x,IV_lookup)
        #println("P3")
        if nest_in_alt_ancestry[a,2] && !flag &&nest_position_in_alt_ancestry[a,1]<nest_position_in_alt_ancestry[a,2]
            #println("P4")
            grad+=Probability(IV_lookup,adjusted_ancestry,i,a+dim1,mu)
        end
    end
    return grad
end

function structural_gradient(b,mu,x,S,A,X,Y)
    dim1,dim2=size(x)
    grad=zeros(dim1,dim2);

    sorted_nest_list,nodes=depth_sorted_nest_list(x)
    #println("nodes",nodes)
    IV_lookup=construct_IV_lookup_table(sorted_nest_list,S,X,A,x,b,mu)

    for root=1:1
        for nest=2:dim1
            grad[root,nest]=gradient_root_to_nest(nest,mu,Y,nodes,IV_lookup)
        end
        for alternative=dim1+1:dim2
            grad[root,alternative]=gradient_root_to_alternative(alternative,mu,A,Y,IV_lookup)
        end
    end
    #println("done")

    for nest_origin=2:dim1
        for nest_destination=2:dim1
            if nest_origin!=nest_destination
                #println("test")
                #println("nest: origin: ",nest_origin, "nest_destination: ", nest_destination)
                grad[nest_origin,nest_destination]=gradient_nest_to_nest(nest_origin,nest_destination,mu,Y,nodes,x,IV_lookup)
            end
        end
    end


    for nest=2:dim1
        for alternative=1+dim1:num_alternatives+dim1
            grad[nest,alternative]=gradient_nest_to_alternative(nest,alternative,mu,A,Y,nodes,x,IV_lookup)
        end
    end

    return -1*grad
end


function full_gradient(z,x,S,A,X,Y)
    dim1,dim2=size(x)
    index_shift=dim1*dim2;
    m,p=size(S)

    grad=zeros(p+dim1+index_shift)
    sorted_nest_list,nodes=depth_sorted_nest_list(x)

    function AD_likelihood(y)
        return -1*likelihood(S,X,A,Y,x,sorted_nest_list,nodes,y[1:p],y[(p+1):end])
    end

    Jacobian_likelihood=z->ForwardDiff.gradient(AD_likelihood,z)
    grad[1:(end-index_shift)]=Jacobian_likelihood(z)
    #grad[(end-index_shift+1):end]=full_gradient_AD(z[1:p],z[(p+1):end],x,S,A,X,Y)[:]
    grad[(end-index_shift+1):end]=structural_gradient(z[1:p],z[(p+1):end],x,S,A,X,Y)[:]
    return grad
end


function f_x(y)
    dim1=(1+(m-2));
    dim2=(1+(m-2)+m);
    index_shift=dim1*dim2;
	return -1*likelihood_x(y[1:2*m],y[(2*m+1):(end-index_shift)],reshape(y[(end-index_shift+1):end],dim1,dim2))
end

################################## Optimization Related functions###################

function solve_NL_subproblem_interval(S,X,A,Y,x,z_t)
    dim1,dim2=size(x)
    m,p=size(S)
    sorted_nest_list,nodes=depth_sorted_nest_list(x)

    function AD_likelihood(y)
        return neg_likelihood(S,X_train,A_train,Y_train,x,y[1:p],y[(p+1):p+dim1])
    end

    global_min, minimisers = minimise( y ->AD_likelihood(y),IntervalBox(0..1,1)× IntervalBox(-1..1,p-1) ×IntervalBox(1..1,1)× IntervalBox(1..3,dim1-1), tol=1e-5 )

    return global_min, minimisers
end



function solve_NL_subproblem_bb(S,X,A,Y,x)
    dim1,dim2=size(x)
    m,p=size(S)

    fixed_index=[1, p+1]

    sorted_nest_list,nodes=depth_sorted_nest_list(x)
    println(nodes)

    function AD_likelihood(y)
        return neg_likelihood(S,X,A,Y,x,sorted_nest_list,nodes,y[1:p],y[(p+1):p+dim1])
    end

    Range_list=[(0.0,0.0)]
    for i=2:p
        push!(Range_list,(-3.0,3.0))
    end

    for i=1:dim1
        if nodes[i]>0
            push!(Range_list,(1.0,10))
        else
            push!(Range_list,(1.0,1.0))
        end
    end


    res= bboptimize(AD_likelihood; SearchRange = Range_list, MaxSteps=40000, Method=:adaptive_de_rand_1_bin_radiuslimited)

    Hessian_likelihood=z->ForwardDiff.hessian(AD_likelihood,z)

    z=best_candidate(res)
    f=best_fitness(res)


    valid_idx=.![idx in fixed_index for idx in collect(1:p+dim1)]
    full_hessian=Hessian_likelihood(z)
    #gen_inv_hessian=Hermitian(pinv(full_hessian[valid_idx,valid_idx]))
    #chol_inv_hessian=cholesky(gen_inv_hessian).L
    #pseudo_var=chol_inv_hessian'*chol_inv_hessian
    std_errs=sqrt.((diag(pinv(full_hessian[valid_idx,valid_idx]))))
    #std_errs=sqrt.((diag(inv(full_hessian[valid_idx,valid_idx]))))

    full_std_errs=zeros(p+dim1)
    full_std_errs[valid_idx]=std_errs


    return z,f,full_std_errs
end

function solve_NL_subproblem(S,X,A,Y,x,z_t;depth_nest_constraint_flag=false)

    dim1,dim2=size(x)
    m,p=size(S)

    fixed_index=[1, p+1]
    sorted_nest_list,nodes=depth_sorted_nest_list(x)

    function JuMP_likelihood(y...)
    	return -1*likelihood(S,X,A,Y,x,sorted_nest_list,nodes,y[1:p],y[(p+1):end])
    end

    function AD_likelihood(y)
    	return -1*likelihood(S,X,A,Y,x,sorted_nest_list,nodes,y[1:p],y[(p+1):p+m-1])
    end
    Hessian_likelihood=z->ForwardDiff.hessian(AD_likelihood,z)
    #Jacobian_likelihood=z->ForwardDiff.gradient(AD_likelihood,z)

    #println("Initializing NL subproblem...")
    #model=Model(solver=IpoptSolver())

    #model=Model(with_optimizer(EAGO.Optimizer))
    #outlev=0
    #solver=NLoptSolver(algorithm=:LD_SLSQP)
    #model=Model(solver=NLoptSolver(algorithm=:LD_SLSQP))
    model=Model(with_optimizer(Ipopt.Optimizer,max_iter=100,tol=0.01,print_level=0))
    #model=Model(with_optimizer(KNITRO.Optimizer,bar_murule=2,algorithm=1,hessopt=2,par_numthreads=8,outlev=0))
    JuMP.register(model, :JuMP_likelihood, p+dim1,JuMP_likelihood, autodiff=true)

    yUpper=ones(p+dim1)
    yLower=zeros(p+dim1)

    for j ∈ 2:p
        yUpper[j]=10
        yLower[j]=-10
    end

    for j∈ 2:dim1
        if nodes[j]>0
            #println("nest ",j,"parent ",nodes[j])
            yUpper[j+p]=10
            yLower[j+p]=1
        else
            yUpper[j+p]=1
            yLower[j+p]=1
            append!(fixed_index,j+p)
        end
    end
    #println("fixed_index:", fixed_index)
    #println("y_lower:", yLower)
    #println("y_upper:", yUpper)


    @variable(model, yLower[j]<=y[j=1:p+dim1]<=yUpper[j])

    flag=true
    for nest in sorted_nest_list
        if flag
            @constraint(model,y[p+nest]<=mu_upper)
            flag=false
        end
        if nodes[nest]>0
            #println(nest," > ",nodes[nest])
            @constraint(model,y[p+nest]>=delta+y[p+nodes[nest]])
            #@constraint(model,y[p+nest]>=y[p+nodes[nest]])
        end
    end

    if depth_nest_constraint_flag==true
        println("Depth constrained nest scales triggered")
        depth_map=construct_depth_map(x)
        max_depth,~=findmax(depth_map)

        for d =1:max_depth
            nests_at_depth=[y for y ∈1:length(depth_map) if depth_map[y]==d]
            for nest in nests_at_depth
                num_nests_at_depth=length(nests_at_depth)
                if num_nests_at_depth>=2
                    for i=2:num_nests_at_depth
                        @constraint(model,y[p+nests_at_depth[i]]==y[p+nests_at_depth[1]])
                    end
                end
            end

        end
    end

    #for j=2:dim1
        #@constraint(model,y[p+j]>=mu_norm)
        #@constraint(model,y[p+j]<=mu_upper)
    #end

    @constraint(model,y[p+1]==mu_norm)
    @constraint(model,y[1]==0)

    #for i=1:p+dim1
        #setvalue(y[i], start_value[i])
    #end


    #objexpr = Expr(:call, :JuMP_likelihood, y...)
    @NLobjective(model, Min, JuMP_likelihood(y...))
    #JuMP.setNLobjective(model, :Min, objexpr)

    #@NLobjective(model, Min, JuMP_likelihood(y...))
    println("Estimating model parameters for tree:")
    print_edges(x)

    optimize!(model)
    #optimize!(model)

    println("Estimation Completed")
    println("-------------------")

    nodes=construct_ancestry(x);
    z=JuMP.value.(y)

    #println(z)

    valid_idx=.![idx in fixed_index for idx in collect(1:p+dim1)]
    full_hessian=Hessian_likelihood(z)
    #gen_inv_hessian=Hermitian(pinv(full_hessian[valid_idx,valid_idx]))
    #chol_inv_hessian=cholesky(gen_inv_hessian).L
    #pseudo_var=chol_inv_hessian'*chol_inv_hessian
    std_errs=zeros(p+dim1-length(fixed_index))
    try
        std_errs=sqrt.(abs.(diag(pinv(full_hessian[valid_idx,valid_idx]))))
    catch
        println("Warning: Hessian not invertible")
    end
    #std_errs=sqrt.((diag(inv(full_hessian[valid_idx,valid_idx]))))

    full_std_errs=zeros(p+dim1)
    full_std_errs[valid_idx]=std_errs



    if length(z_t)>0
        b_t=z_t[1:p]
        mu_t=z_t[p+1:end]

        println("------Results------")
        println("ASCs")
        println("___________________")
        println("|-----|","True  ","|","Est. |","std. |","t-stat. |")
        for i=1:m
            println("|ASC_",i,"|",b_t[i],"   |",round(z[i],digits=2),"  |",round(full_std_errs[i],digits=2),"  |",round(z[i]/full_std_errs[i],digits=2))
        end
        println("-------------------")
        println("Betas")
        println("___________________")
        println("|---|","True  ","|","Est. |","std. |","t-stat. |")

        for i=m+1:p
            println("|β_",i-m,"|",b_t[i],"  |",round(z[i],digits=2),"  |",round(full_std_errs[i],digits=2),"  |",round(z[i]/full_std_errs[i],digits=2))
        end

        println("-------------------")
        println("Nest Scales")
        println("___________________")
        println("|---|","True ","|","Est.  |","std.  |","t-stat.|")
        for i=1:dim1
            if nodes[i]>0
                println("|μ_",i,"|",mu_t[i],"  |",round(z[p+i],digits=1),"   |",round(full_std_errs[p+i],digits=2),"  |",round((z[p+i]-mu_t[nodes[i]])/full_std_errs[p+i],digits=2))
            end
        end
        println("-------------------")
        println("Log-likelihood")
        println("___________________")
        println("|---|","True ","|","Est.|")
        println("|LL| ",round(AD_likelihood(z_t),digits=2),"  |",round(AD_likelihood(z),digits=2),"|")
    end
    return z, full_std_errs, getobjectivevalue(model)
end

function construct_sub_mnl_inputs(R,S,X,A,C,children)
    #println(children)
    n,~=size(A)
    S_mnl=S[children,:]
    R_mnl=R[children,:]
    relevant_paramters=(sum(S_mnl,dims=1)'.>0)[:,1]
    S_mnl=S_mnl[:,relevant_paramters]
    R_mnl=R_mnl[:,relevant_paramters]

    #w,v=size(C)
    #print(w)
    #print(v)

    C_mnl=C[:,children]
    relevant_observations=(sum(C_mnl,dims=2).>0)[:,1]
    C_mnl=C_mnl[relevant_observations,:]

    #relevant_observations=[Y[i] ∈ children for i=1:n]
    #Y_mnl=Y[relevant_observations]

    A_mnl=A[relevant_observations,children]
    X_mnl=X[relevant_observations,relevant_paramters]


    return R_mnl,S_mnl,X_mnl,A_mnl,C_mnl,relevant_paramters
end

function prune_tree(nest,children,index_shift,γ,S,X,A,C)
    println("pruning nest: ",nest)
    m,p=size(S)
    n,m=size(A)


    num_children=length(children)

    logsums=zeros(n)
    for i=1:n
        for j=1:m
            if A[i,j]>0
                logsums[i]+=exp(sum(γ[k]*S[j,k]*X[i,k] for k=1:p))
            end
        end
        logsums[i]=log(logsums[i])

        A[i,nest]=any(A[i,children].>0)*1
        A[i,children]=zeros(num_children)

        C[i,nest]=any(C[i,children].>0)*1
        C[i,children]=zeros(num_children)


    end

    X[:,index_shift+nest]=logsums

    return S,X,A,C
end

function estimate_sub_mnl(R,S,X,A,C)
    n,m=size(A)
    ~,p=size(S)

    println(p)
    println(n)
    println(m)

    function likelihood_mnl(y)
        ll=Array{Any}(undef,n)
        norm_util=Array{Any}(undef,n,m)
        Threads.@threads for i=1:n
            for j=1:m
                if A[i,j]>0
                    norm_util[i,j]=exp(sum(y[k]*S[j,k]*X[i,k] for k=1:p))
                else
                    norm_util[i,j]=0
                end
            end
            for j=1:m
                if C[i,j]>0
                    ll[i]=log(norm_util[i,j])-log(epsilon+sum(norm_util[i,:]))
                end
            end
            #println(ll)
            #ll[i]=sum(y[k]*S[Y[i],k]*X[i,k] for k=1:p)-
            #println(ll[i])
        end
        #print(ll)
        return -1*sum(ll)
    end

    function JuMP_likelihood_mnl(y...)
        return likelihood_mnl(y[1:p])
    end

    if true
        #model=Model(solver=IpoptSolver())
        model=Model(solver=KnitroSolver(bar_murule=2, convex=1,algorithm=2,par_numthreads=8,maxit=5000))
        JuMP.register(model, :JuMP_likelihood_mnl, p,JuMP_likelihood_mnl, autodiff=true)

        @variable(model, y[1:p])
        #@constraint(model,y[1]==0)

        for i=1:p
            if sum(R[:,i])>0
                println("triggering restriction constraints")
                @constraint(model,y[i]<=1)
                @constraint(model,y[i]>=0)
            end
        end

        #objexpr = Expr(:call, :JuMP_likelihood_mnl, y...)
        #JuMP.setNLobjective(model, :Min, objexpr)
        @NLobjective(model, Min, JuMP_likelihood_mnl(y...))
        optimize!(model)

        return JuMP.value.(y)
    else
        lower=zeros(p)
        upper=zeros(p)
        initial=zeros(p)
        for i=2:p
            if sum(R[:,i])>0
                upper[i]=1
                initial[i]=0.5
            else
                lower[i]=-10
                upper[i]=10
            end
        end
        #func = TwiceDifferentiable(y -> likelihood_mnl(y),
                         #  ones(p); autodiff=:forward);

        inner_optimizer =GradientDescent()


        opt = Optim.optimize(likelihood_mnl, lower, upper, initial, Fminbox(inner_optimizer); autodiff = :forward)

        #Hessian_likelihood=z->ForwardDiff.hessian(AD_likelihood,z

        return Optim.minimizer(opt)
    end
end

function sequential_nested_logit(S,X,A,C,x)
    m,p=size(S)
    n,~=size(A)
    dim1,dim2=size(x)

    μ=ones(m-1)
    β=zeros(p)

    index_shift=p+dim1


    S_expanded=zeros(dim2,p+2*dim1)
    S_expanded[1:dim2,1:dim2]=Matrix{Float64}(I, dim2, dim2)
    S_expanded[1:dim1,dim1+p+1:end]=Matrix{Float64}(I, dim1, dim1)
    S_expanded[dim1+1:end,dim2+1:dim2+p-m]=S[:,m+1:end]

    R=zeros(dim2,p+2*dim1)
    R[1:dim1,dim1+p+1:end]=Matrix{Float64}(I, dim1, dim1)

    x_expanded=copy(x)

    X_expanded=hcat(ones(n,dim1),X,zeros(n,dim1))
    A_expanded=hcat(zeros(n,dim1),A)
    C_expanded=hcat(zeros(n,dim1),C)

    γ=zeros(p+2*dim1)

    sorted_nest_list,nodes=depth_sorted_nest_list(x)
    depth_map=construct_depth_map(x)

    #println(nodes)
    #print_edges(x)

    for nest ∈ sorted_nest_list
        println("Estimating: Nest# ",nest," @ depth ",depth_map[nest])
        children=get_children(nest,x)
        #println(children)
        R_mnl,S_mnl,X_mnl,A_mnl,C_mnl,relevant_parameters=construct_sub_mnl_inputs(R,S_expanded,X_expanded,A_expanded,C_expanded,children)
        γ[relevant_parameters]=estimate_sub_mnl(R_mnl,S_mnl,X_mnl,A_mnl,C_mnl)
        S_expanded,X_expanded,A_expanded,C_expanded=prune_tree(nest,children,index_shift,γ,S_expanded,X_expanded,A_expanded,C_expanded)
    end

    for nest ∈ sorted_nest_list[end-1:-1:1]
        parent=nodes[nest]
        μ[nest]=μ[parent]/γ[index_shift+nest]
        γ[nest]=γ[nest]/μ[parent]
    end

    for nest ∈ 1:dim1
        for child_nest ∈ nest:dim1
            if x[nest,child_nest]>0
                γ[child_nest]+=γ[nest]
            end
        end
    end



    for alternative ∈ dim1+1:dim2
        parent=nodes[alternative]
        relevant_parameters=S_expanded[alternative,:].>0
        β[S[alternative-dim1,:].>0]=γ[relevant_parameters]/μ[parent]
        β[alternative-dim1]+=γ[parent]
    end

    for alternative ∈ m:-1:1
        β[alternative]+=-1*β[1]
    end

    return β,μ

end

function solveMasterProblem(hyperplanes,upperbounds_list,visited_x,visited_solutions,prev_opt_solution,visited_sol_count, max_iter,m,p,depth,nests)
    #println("Solving Master Problem...")
    dim1=1+(m-2)
    dim2=dim1+m

    #println(visited_sol_count)

    MasterProblem = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LazyConstraints" => 1))

    @variable(MasterProblem, x[1:dim1,1:dim2], Bin)
    @variable(MasterProblem, μ[1:dim1]>=1)
    @variable(MasterProblem, 0<=y[1:m-2]<=1,Bin)
    @variable(MasterProblem, β[1:p])
    @variable(MasterProblem, η)

    #@variable(MasterProblem,>=0)
    @variable(MasterProblem, δ, Bin)

    ϵ=1
    #bigM=10^4

    λ=ϵ*exp(log(1/ϵ)*visited_sol_count/max_iter)
    #λ=0

    #@objective(MasterProblem , Min, η+λ*ξ)
    #@objective(MasterProblem , Min, η)

    if visited_sol_count==0
        @objective(MasterProblem , Min, η)
        @constraint(MasterProblem,η>=0)
    else
        #println("Trigger")
        @variable(MasterProblem,ξ[1:visited_sol_count]>=0)
        @objective(MasterProblem , Min, η+λ*sum(ξ[j] for j=1:visited_sol_count))
        for i=1:visited_sol_count
            Δ=vcat(β,μ,x[:])-vcat(visited_solutions[:,i],visited_x[:,:,i][:])
            @constraint(MasterProblem, η>=-ξ[i]+dot(hyperplanes[:,i],Δ)+upperbounds_list[i])
            #@constraint(MasterProblem, η>=dot(hyperplanes[:,i],Δ)+upperbounds_list[i])

            @constraint(MasterProblem, η<=upperbounds_list[i])

            #cut off Previously visited solutions
            @constraint(MasterProblem,sum(sum((2*visited_x[j,k,i]-1)*x[j,k] for j=1:dim1) for k=1:dim2) <= sum(sum((visited_x[j,k,i]) for j=1:dim1) for k=1:dim2)-1)
        end
    end

    #@constraint(MasterProblem,η>=0)

    if length(prev_opt_solution)>0
        #println(prev_opt_solution)
        depth_map=construct_depth_map(prev_opt_solution)
        nodes=construct_ancestry(prev_opt_solution)
        max_depth,~=findmax(depth_map)

        for d =0:max_depth
            #println("Depth: ",d)
            nests_at_depth=[y for y ∈1:length(depth_map) if depth_map[y]==d]
            for nest in nests_at_depth
                #println("Nest: ",nest)
                parent=nodes[nest]
                if parent>0
                    @constraint(MasterProblem,x[parent,nest]==1)
                end
                children=get_children(nest,prev_opt_solution)
                for child in children
                    #println("Child: ",child)
                    if d==max_depth
                        #println("flag1")
                        for other_nests in setdiff(1:dim1,nest)
                            if nodes[other_nests]>0 || other_nests==1
                                @constraint(MasterProblem, x[other_nests,child]==0)
                            end
                        end
                        #@constraint(MasterProblem, x[nest,child]+sum(x[nest,k]*x[k,child] for k=2:dim1)==1)
                    else
                        #println("flag2")
                        @constraint(MasterProblem, x[nest,child]==1)
                    end

                end
            end
        end

    end



    #@constraint(MasterProblem,x[1,8]==1)
    #@constraint(MasterProblem,sum(x[1,j] for j=1:dim2)==2)


    # Constraint total number of nests
    @constraint(MasterProblem,sum(y[i] for i=1:dim1-1)==nests)

    #Constraint Max depth
    if length(prev_opt_solution)==0
        longest_branch = zero(AffExpr)
        for i=1:dim1-1
            #global longest_branch
            longest_branch+=x[i,i+1]
        end
        @constraint(MasterProblem,longest_branch==depth)
    end

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
        @constraint(MasterProblem,sum(x[j,i+1] for j=1:dim1)==y[i])
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

    for i=1:p
        @constraint(MasterProblem,β[i]<=10)
        @constraint(MasterProblem,β[i]>=-10)
    end
    @constraint(MasterProblem,β[1]==0)
    @constraint(MasterProblem,μ[1]==mu_norm)

    #Graph is a forest
    @constraint(MasterProblem,sum(sum(x[i,j] for i=1:dim1) for j=1:dim2)==(sum(y[i] for i=1:m-2)+m+1)-1)




    function cycle_detection(cb_data,x_val)
        #println("----\nInside cycle detection callback")

        A=convert(Array{Int64},vcat(x_val,zeros(m,dim2)))
        Tree=DiGraph(A.|A')

        if is_cyclic(Tree)
            cycles= simplecycles(Tree)
            for cycle in cycles
                if length(cycle)<3
                    continue
                else
                    arcs_in_cycle = zero(AffExpr)
                    nodes_in_cycle= zero(AffExpr)
                    cycle_length=length(cycle)
                    #println("A cycle of length ",cycle_length," has been detected")
                    #println("cycle: ",cycle)
                    #println("A",A)
                    #println("x",x_val)
                    #println("x",JuMP.value.(x))

                    nodes_in_cycle = zero(AffExpr)


                    for i=1:cycle_length
                        if 1<cycle[i]<=dim1
                            nodes_in_cycle+=y[cycle[i]-1]
                        else
                            nodes_in_cycle+=1
                        end
                    end

                    arcs_in_cycle = zero(AffExpr)

                    if cycle[1]<=dim1
                        arcs_in_cycle+=x_val[cycle[1],cycle[2]]*x[cycle[1],cycle[2]]
                    end

                    if cycle[2]<=dim1
                        arcs_in_cycle+=x_val[cycle[2],cycle[1]]*x[cycle[2],cycle[1]]
                    end


                    #arcs_in_cycle+=x_val[cycle[1],cycle[2]]*x[cycle[1],cycle[2]]+x_val[cycle[2],cycle[1]]*x[cycle[2],cycle[1]]
                    for i=2:cycle_length-1
                        if cycle[i]<=dim1
                            arcs_in_cycle+=x_val[cycle[i],cycle[i+1]]*x[cycle[i],cycle[i+1]]
                        end

                        if cycle[i+1]<=dim1
                            arcs_in_cycle+=x_val[cycle[i+1],cycle[i]]*x[cycle[i+1],cycle[i]]
                        end

                    end

                    if cycle[cycle_length]<=dim1
                        arcs_in_cycle+=x_val[cycle[cycle_length],cycle[1]]*x[cycle[cycle_length],cycle[1]]
                    end

                    if cycle[1]<=dim1
                        arcs_in_cycle+=x_val[cycle[1],cycle[cycle_length]]*x[cycle[1],cycle[cycle_length]]
                    end
                    #println("Adding cycle elimination cut")
                    con=@build_constraint( arcs_in_cycle <= nodes_in_cycle-1)
                    MOI.submit(MasterProblem, MOI.LazyConstraint(cb_data), con)

                end

            end
        else
            #println("No cycles detected!")
            return
        end

    end

    function equivalence_cut(cb_data,x_val)
        #println("----\nInside Equivalence cut generator")

        flag=false
        for i=1:visited_sol_count
            flag =determine_equivalence(x_val,visited_x[:,:,i])
            if flag==true
                con=@build_constraint(sum(sum((2*x_val[j,k]-1)*x[j,k] for j=1:dim1) for k=1:dim2) <= sum(sum((x_val[j,k]) for j=1:dim1) for k=1:dim2)-1)
                MOI.submit(MasterProblem, MOI.LazyConstraint(cb_data), con)
            end
        end
        #println("No more cuts!")
        return
    end

    function depth_cut(cb_data,x_val)
        #println("Inside pruning callback")

        #current_x=round.(callback_value(cb_data,x),digits=1)

        nodes=construct_ancestry(x_val)
        for j=dim1+1:dim2
            #println("looking at alternative ",j)
            path=path_from_root(j,nodes)
            #println("path computed")
            l=length(path)
            depth_path=l-2
            if depth_path>depth
                #println("Maximum depth Exceeded")
                branch = zero(AffExpr)
                for i=1:l-2
                    branch+=x[path[i],path[i+1]]
                end
                #println("Adding prunning constraint")
                con=@build_constraint(branch<=depth)
                MOI.submit(MasterProblem, MOI.LazyConstraint(cb_data), con)
            end
        end
        #println("Tree pruned!")
        return
    end

    function lazy_cut_wrapper(cb_data)
        #add constraints at integer nodes
        x_val=zeros(dim1,dim2)
        for i ∈1:dim1
            for j∈1:dim2
                x_val[i,j]=callback_value(cb_data,x[i,j])
            end
        end
        if is_valid_tree(x_val)
            equivalence_cut(cb_data,x_val)
            depth_cut(cb_data,x_val)
            cycle_detection(cb_data,x_val)
        end
    end


    #MOI.set(MasterProblem, MOI.LazyConstraintCallback(), equivalence_cut)
    #MOI.set(MasterProblem, MOI.LazyConstraintCallback(), depth_cut)
    MOI.set(MasterProblem, MOI.LazyConstraintCallback(), lazy_cut_wrapper)
    #addlazycallback(MasterProblem,cycle_detection)
    #addlazycallback(MasterProblem,equivalence_cut)
    #addlazycallback(MasterProblem,depth_cut)
    optimize!(MasterProblem)
    status=termination_status(MasterProblem)

    #println("ξ ",JuMP.value.(ξ))
    #println("λ ",λ)
    #println("μ",JuMP.value.(μ))


    if status== MOI.OPTIMAL
        x_sol=JuMP.value.(x)
        y_sol=JuMP.value.(y)

        for j=2:dim1
            x_sol[j,:]=y_sol[j-1]*x_sol[j,:]
        end

        x_val=(round.(x_sol,digits=1).>0)*1
        η=JuMP.value.(η)
    else
        η=1e9
        x_val=zeros(dim1,dim2)
    end



    #x_val=convert(Array{Int64},x_val)


    #println("solution: ", x_val)
    #println("solution y: ", JuMP.value.(y))

    #round.(JuMP.value.(x),digits=1)
    #println(x_val)
    #println(is_valid_tree(x_val))

    return status, x_val, η
end


function optimal_logit(S,X,A,Y,X_val,A_val,Y_val,z_t,prev_opt_solution,depth,nests,max_iter)
    count=0
    m,p=size(S)
    dim1=1+(m-2)
    dim2=dim1+m

    upperbounds_list=(1e9)*ones(max_iter)
    lowerbounds_list=zeros(max_iter)
    visited_x=zeros(dim1,dim2,max_iter)
    hyperplanes=zeros(p+dim1+dim1*dim2,max_iter)
    visited_solutions=zeros(p+dim1,max_iter)
    standard_errors=zeros(p+dim1,max_iter)
    validation_list=(1e9)*ones(max_iter)


    for i=1:max_iter
        status, visited_x_temp,lowerbounds_temp=solveMasterProblem(hyperplanes,
        upperbounds_list,visited_x,visited_solutions,prev_opt_solution,count,max_iter,m,p,depth,nests)

        if status!= MOI.OPTIMAL
            print("Infeasibility")
            break
        end
        visited_x[:,:,i]=visited_x_temp
        lowerbounds_list[i]=lowerbounds_temp
        visited_solutions[:,i], standard_errors[:,i], upperbounds_list[i]=solve_NL_subproblem(S,X,A,Y,visited_x[:,:,i],[],depth_nest_constraint_flag=false)
        sorted_nest_list,nodes=depth_sorted_nest_list(visited_x[:,:,i])
        validation_list[i]=neg_likelihood(S,X_val,A_val,Y_val,visited_x[:,:,i],sorted_nest_list,nodes,visited_solutions[1:p,i],visited_solutions[p+1:end,i])
        #~,~,validation_list[i]=solve_NL_subproblem(S,X_val,A_val,Y_val,visited_x[:,:,i],[])
        hyperplanes[:,i]=full_gradient(visited_solutions[:,i],visited_x[:,:,i],S,A,X,Y)
        count+=1
    end
    return validation_list,upperbounds_list,lowerbounds_list,visited_x,visited_solutions,standard_errors,hyperplanes,count
end

function nest_tree_search(S,X,A,Y,X_val,A_val,Y_val,max_iter)

    m,~=size(S)
    dim1=m-1
    dim2=dim1+m

    optimal_trees=zeros(dim1,dim2,m-2,m-2)
    optimal_ll=(1e9)*ones(m-2,m-2)
    #validation_ll=zeros(m-2)
    optimal_params=zeros(p+dim1,m-2,m-2)
    optimal_std_errs=zeros(p+dim1,m-2,m-2)

    counts=zeros(m-2,m-2)

    prev_opt_solution=[]
    f_depth_best=1e9
    for depth=1:m-2
        f_nest_best=1e9
        for nests=depth:m-2
            println("Nests: ", nests," @ depth: ", depth)
            validation_list,upperbounds_list,lowerbounds_list,visited_x,visited_solutions,
            standard_errors,hyperplanes,counter=optimal_logit(S,X,A,Y,X_val,A_val,Y_val,[],prev_opt_solution,depth,nests,max_iter)

            #if counter==0
                #break
            #end

            counts[depth,nests]=counter
            #optimal_ll[depth,nests],best_idx=findmin(validation_list)
            optimal_ll[depth,nests],best_idx=findmin(upperbounds_list)
            if optimal_ll[depth,nests]<f_nest_best
                f_nest_best=optimal_ll[depth,nests]
            else
                break
            end

            #optimal_ll[depth],best_idx=findmin(upperbounds_list[upperbounds_list.>0])
            optimal_trees[:,:,depth,nests]=visited_x[:,:,best_idx]
            optimal_params[:,depth,nests]=visited_solutions[:,best_idx]
            optimal_std_errs[:,depth,nests]=standard_errors[:,best_idx]
        end
        #if sum(optimal_ll[depth,:].>0)*1==0
            #break
        #end
        f_min,best_of_best_idx=findmin(optimal_ll[depth,:])
        prev_opt_solution=optimal_trees[:,:,depth,best_of_best_idx]
        if f_min<f_depth_best
            f_depth_best=f_min
        else
            break
        end

        if f_min==1e9
            break
        end
        println("Optimal tree at current depth:")
        println("************************************")
        println("************************************")
        println("************************************")
        print_edges(prev_opt_solution)
        println("************************************")
        println("************************************")
        println("************************************")
        draw_tree(prev_opt_solution,optimal_params[p+1:end,depth,best_of_best_idx])
    end


    return optimal_trees,optimal_ll,optimal_params,optimal_std_errs,counts
end

function solve_NLSLP(S,X,A,Y,X_val,A_val,Y_val,max_iter)
    m,p=size(S)
    dim1=m-1
    dim2=dim1+m

    optimal_trees=zeros(2,dim1,dim2,m-2,m-2)
    optimal_ll=(1e9)*ones(2,m-2,m-2)
    optimal_params=zeros(2,p+dim1,m-2,m-2)
    optimal_std_errs=zeros(2,p+dim1,m-2,m-2)

    counts=zeros(m-2,m-2)
    for depth=1:m-2
        for nests=depth:1:m-2
            println("Nests: ", nests," @ depth: ", depth)
            validation_list,upperbounds_list,lowerbounds_list,visited_x,visited_solutions,
            standard_errors,hyperplanes,counter=optimal_logit(S,X,A,Y,X_val,A_val,Y_val,[],[],depth,nests,max_iter)


            counts[depth,nests]=counter
            optimal_ll[1,depth,nests],best_idx=findmin(upperbounds_list)
            optimal_trees[1,:,:,depth,nests]=visited_x[:,:,best_idx]
            optimal_params[1,:,depth,nests]=visited_solutions[:,best_idx]
            optimal_std_errs[1,:,depth,nests]=standard_errors[:,best_idx]

            optimal_ll[2,depth,nests],best_idx=findmin(validation_list)
            optimal_trees[2,:,:,depth,nests]=visited_x[:,:,best_idx]
            optimal_params[2,:,depth,nests]=visited_solutions[:,best_idx]
            optimal_std_errs[2,:,depth,nests]=standard_errors[:,best_idx]

        end
    end

    best_val_ll,best_val_idx=findmin(optimal_ll[2,:,:])
    best_train_ll,best_train_idx=findmin(optimal_ll[1,:,:])
    opt_val_tree=optimal_trees[2,:,:,best_val_idx]
    opt_val_params=optimal_params[2,:,best_val_idx]
    opt_val_std_errs=optimal_std_errs[2,:,best_val_idx]

    opt_train_tree=optimal_trees[1,:,:,best_train_idx]
    opt_train_params=optimal_params[1,:,best_train_idx]
    opt_train_std_errs=optimal_std_errs[1,:,best_train_idx]

    println("Optimal VALIDATION tree (validation loglikekihood= ",best_val_ll,"): ")
    println("************************************")
    println("************************************")
    print_edges(opt_val_tree)
    draw_tree(opt_val_tree,opt_val_params[p+1:end])
    println("************************************")
    println("************************************")
    println("************************************")
    println("Model Parameters and (standard errors)")
    for i=1:p
        println("|β_",i,"|",opt_val_params[i],"  (",round(opt_val_std_errs[i],digits=2),")")
    end
    println("Nest scale parameters and (standard errors)")
    nodes=construct_ancestry(opt_val_tree)
    for i=1:dim1
        if nodes[i]>0
            println("|μ_",i,"|",opt_val_params[i+p],"  (",round(opt_val_std_errs[p+i],digits=1),")")
        end
    end


    println("Optimal TRAINING tree (training loglikekihood= ",best_train_ll,"): ")
    println("************************************")
    println("************************************")
    print_edges(opt_train_tree)
    draw_tree(opt_train_tree,opt_train_params[p+1:end])
    println("************************************")
    println("************************************")
    println("************************************")
    println("Model Parameters and (standard errors)")
    for i=1:p
        println("|β_",i,"|",opt_train_params[i],"  (",round(opt_train_std_errs[i],digits=2),")")
    end
    println("Nest scale parameters and (standard errors)")
    nodes=construct_ancestry(opt_train_tree)
    for i=1:dim1
        if nodes[i]>0
            println("|μ_",i-1,"|",opt_train_params[i+p],"  (",round(opt_train_std_errs[p+i],digits=1),")")
        end
    end


    for l ∈ 1:2
        for i ∈ 1:m-2
            for j ∈ 1:m-2
                if optimal_ll[l,i,j]==1e9
                    optimal_ll[l,i,j]=NaN
                end
            end
        end
    end

    return optimal_trees,optimal_ll,optimal_params,optimal_std_errs,counts

end

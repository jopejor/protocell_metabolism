using Distributions
using Combinatorics
using Iterators
using ClusterManagers
using DataFrames

#Functions 
function gmean(A::Array{Int64,1})
	# geometric average of a vector
	value=prod(A)^(1/length(A))
	return value
end

function xindicator(s1,s2)

	if s1 == 0 
		value=0
	elseif s2 ==0
		value=0
	elseif s2 > s1
		value=s1/s2
	elseif s1 >= s2
		value=1
	end		
return value
end

function Metabolism(x)
        flux=x[1]+1.5*x[2]*2*xindicator(x[1],x[2]*2)+0.00000001
return flux
end



function states(state_vec)

	state_zero = collect(permutations(state_vec[1]))
	for i=1:length(state_vec)
		additional=collect(permutations(state_vec[i]))
		state_zero=append!(state_zero,additional)
	end

	state_zero=unique(state_zero)

	return state_zero
end

function multinomial2(k)
    s = 0
    result = 1
    @inbounds for i in k
        s += i
        result *= binomial(BigInt(s), i)
    end
    result
end


function development(state_zero,N,Vec,mu)
	q=(1-mu)^100
	p=deepcopy(state_zero)/N
 	NewVec=zeros(length(state_zero))
    for j=1:length(state_zero)
	    for i=1:length(state_zero)
	    	prob=multinomial2(state_zero[j])*prod((q.*p[i][1:end-1]).^state_zero[j][1:end-1])*(q*p[i][end]+(1-q))^state_zero[j][end]
	    	NewVec[j]=NewVec[j]+prob*Vec[i]
	    end
	end

	return NewVec
end


function fitness(state_zero,NewVec,Metabolism)

	fit_vec=zeros(length(state_zero))

		for i=1:length(state_zero)
			fit_vec[i]=Metabolism(state_zero[i])
	    end

	vec=fit_vec.*NewVec
	Norm_Vec=vec/sum(vec)
	lamda=sum(vec)

	return Norm_Vec,lamda
end

function coexitence_fraction(state_zero,Vec,seqtypes)

		coex_states=zeros(length(state_zero))

		for i=1:length(state_zero)
			if all(x->x>0,state_zero[i][1:end-1])
			coex_states[i]=1
			else
				coex_states[i]=0
			end
	    end
	fraccio=sum(coex_states.*Vec)

    return fraccio
end


function perron_iter(N,seqtypes,mu)

	tol=0.001
	state_vec=collect(partitions(N+seqtypes,seqtypes))-1
	state_zero=states(state_vec)
	vec0=1/length(state_zero)*ones(length(state_zero))
	difference=1
	Vec=vec0

	while difference > tol
		NewVec=development(state_zero,N,Vec,mu)
		NewVec=fitness(state_zero,NewVec,Metabolism)[1]
		lambda=fitness(state_zero,NewVec,Metabolism)[2]
		difference=sum(abs(NewVec-Vec))
		Vec=NewVec
	end
	fraccio=coexitence_fraction(state_zero,Vec,seqtypes)
	lambda=fitness(state_zero,Vec,Metabolism)[2]


    file_prefix= join([seqtypes,N,mu],"_")


	df4 = DataFrame()
    df4[:lambda]=lambda
    df4[:coex]=fraccio

    distribucions_name=join([file_prefix, ".csv"],"_")
    writetable(distribucions_name, df4)

end








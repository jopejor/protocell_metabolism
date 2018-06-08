using Distributions
using Iterators
using ClusterManagers
using DataFrames

require("types.jl")

#-----------------------------
# Functions For Protcells
#-----------------------------


function Metabolism(x)
    #flux function, substitute with the flux to be studied
   flux=x[1]+2*min(x[1],2*x[2])
return flux
end


function InitialPop(popsize,seqtypes,N,Metabolism)

    Population_zero = Array(Protocell,popsize)
    equipartite=round(N/(seqtypes-1))
    equipartite=convert(Integer,equipartite)
    Founder=zeros(Integer,seqtypes-1)
    fill!(Founder,equipartite)
    push!(Founder,0)

    if sum(Founder)>N
    	Diff=N-sum(Founder)
    	where=rand(1:seqtypes-1)
        Founder[where]=Founder[where]+Diff
    elseif sum(Founder)<N
    	Diff=N-sum(Founder)
    	where=rand(1:seqtypes-1)
        Founder[where]=Founder[where]+Diff
    end
    developmenta_init=0
    initialProt=Protocell(Founder,0.1,Metabolism)
    for i=1:popsize
        Population_zero[i]=deepcopy(initialProt)
    end
    Initial_population=Population(Population_zero)

    return Initial_population
end


function SequenceDynSelf(protocell,mu,L,N)

	q=(1-mu)^L
	total=sum(protocell.composition)
	temps=0

    while total != 2*N

        #Pick the  sequence type

        sec_freq=protocell.composition/total
        d=Categorical(sec_freq)
        R=rand(d)

        bi=Binomial(1,q)
        test=rand(bi)

        if test == 1
            protocell.composition[R]=protocell.composition[R]+1
        elseif test == 0
            protocell.composition[end]=protocell.composition[end]+1
        end
        total=sum(protocell.composition)
    end


    return protocell
end


function ProtocellSplit(protocell,N)

    total=sum(protocell.composition)
    Daughter_1=zeros(Integer,length(protocell.composition))
    Daughter_2=zeros(Integer,length(protocell.composition))
    
    sum_daughter=sum(Daughter_1)
    
    while sum_daughter != N

        total=sum(protocell.composition)
        protocell_freq=protocell.composition/total
        d=Categorical(protocell_freq)
        R=rand(d)
     
        while protocell.composition[R] == 0
            d=Categorical(protocell_freq)
            R=rand(d)
        end

        Daughter_1[R]=Daughter_1[R]+1
        sum_daughter=sum(Daughter_1)

        protocell.composition[R]=protocell.composition[R]-1
    end

    Daughter_1=Daughter_1
    Daughter_2=protocell.composition

    prot_1=Protocell(Daughter_1,protocell.devtime,Metabolism)
    prot_2=Protocell(Daughter_2,protocell.devtime,Metabolism)

    return prot_1,prot_2
end


function Fitness(poblacio)

    fitness_vec=zeros(length(poblacio.individuals))

    for i=1:length(poblacio.individuals)
        if poblacio.individuals[i].devtime != 0
            poblacio.individuals[i].devtime=poblacio.individuals[i].metabolism(poblacio.individuals[i].composition)
            fitness_vec[i]=poblacio.individuals[i].devtime
        else
            fitness_vec[i]=0
        end          
    end
    return fitness_vec
end


#-----------------------------
# Population Dynamics
#-----------------------------

function ForwardGen(poblacio,mu,L,N,popsize)

    NewPop = Array(Protocell,2*length(poblacio.individuals))
 

    #Synchronous replication

    for i=1:length(poblacio.individuals)
       SequenceDynSelf(poblacio.individuals[i],mu,L,N)
    end

    for i=1:length(poblacio.individuals)
        Sons=ProtocellSplit(poblacio.individuals[i],N)
        NewPop[i]=Sons[1]
        NewPop[end-i+1]=Sons[2]
    end

    poblacio=Array(Protocell,popsize)
    poblacio=Population(poblacio)

    #Selective filter
    NewPop=Population(NewPop)

    if length(NewPop.individuals)<=popsize
        poblacio=NewPop
    else 
        for i=1:popsize
            fitness_vec=Fitness(NewPop)
            total=sum(fitness_vec)
            fitness_freq=fitness_vec/total
            d=Categorical(fitness_freq)
            R=rand(d)
            poblacio.individuals[i]=deepcopy(NewPop.individuals[R])
            NewPop.individuals[R].devtime=0
        end
    end
    return poblacio
end

function unwrapper(poblacio) 

    x = Array(Vector,length(poblacio.individuals))

    for i=1:length(poblacio.individuals)
        x[i]=poblacio.individuals[i].composition
    end

    return x
end


function mean_sequence(poblacio,indexgene) 

    x = zeros(0)

    for i=1:length(poblacio.individuals)
        if poblacio.individuals[i].composition[indexgene] >= 0
        push!(x,poblacio.individuals[i].composition[indexgene])
    end
    end
    
    if length(x)== 0
        x=0
    end

    men_gene=mean(x)
    std_gene=std(x)

    return men_gene,std_gene
end

function mean_omega(poblacio) 

    x = zeros(0)

    for i=1:length(poblacio.individuals)
        if poblacio.individuals[i].composition[end] >= 0
        push!(x,poblacio.individuals[i].composition[end])
    end
    end
    
    if length(x)== 0
        x=0
    end

    men_gene=mean(x)
    std_gene=std(x)

    return men_gene,std_gene
end

function coexistence_fraccio(poblacio,seqtypes)

    ihasall=0

        for i=1:length(poblacio.individuals)
            if all(x->x>0,poblacio.individuals[i].composition[1:end-1])
                ihasall=ihasall+1
            end
        end

    fraccio=ihasall/length(poblacio.individuals)

    return fraccio
end


function mean_ploidy(poblacio) 

    x = zeros(0)

    for i=1:length(poblacio.individuals)
        if all(x->x>0,poblacio.individuals[i].composition[1:end-1])
        n_e=sum(poblacio.individuals[i].composition[1:end-1])
        push!(x,n_e)
        end
    end
    
    if length(x)== 0
        x=0
    end

    men_gene=mean(x)
    std_gene=std(x)

    return men_gene,std_gene
end



function SimPop(N,mu)

    #Parameters

    seqtypes = 3
    NumGens = 5000
    popsize=1000
    L=100
    

    
    #generate population with founder

    poblacio=InitialPop(popsize,seqtypes,N,Metabolism)

    #evolve population for NumGens storing mean population fitness

    fraccio_coexistent=zeros(NumGens)
    fitness_promig=zeros(NumGens)
    fitness_sd=zeros(NumGens)
    omega_promig=zeros(NumGens)
    omega_sd=zeros(NumGens)
    effective=zeros(NumGens)
    effective_sd=zeros(NumGens)



    for i=1:NumGens

      poblacio=ForwardGen(poblacio,mu,L,N,popsize)

        if length(poblacio.individuals) == 0
            break
        end 

        fraccio_coexistent[i]=coexistence_fraccio(poblacio,seqtypes)
        fitness_promig[i]=mean(Fitness(poblacio))
        fitness_sd[i]=std(Fitness(poblacio))
        omega_promig[i]=mean_omega(poblacio)[1]
        omega_sd[i]=mean_omega(poblacio)[2]
        effective[i]=mean_ploidy(poblacio)[1]
        effective_sd[i]=mean_ploidy(poblacio)[2]
    end
    
    poblacio_save=unwrapper(poblacio)
    #Saving measurements

    file_prefix= join([mu,seqtypes,N],"_")

    poblacio_name=join([file_prefix, "poblacio.csv"],"_")
    writecsv(poblacio_name,poblacio_save)

    coexistence_name=join([file_prefix, "coexistence_frac.csv"],"_")
    writecsv(coexistence_name,fraccio_coexistent)


    omega_name=join([file_prefix, "omegas.csv"],"_")
    omegadf = DataFrame()
    omegadf[:mean_omega]=omega_promig
    omegadf[:sd_omega]=omega_sd
    writetable(omega_name, omegadf)

    fitness_name=join([file_prefix, "fitness.csv"],"_")
    fitnessdf = DataFrame()
    fitnessdf[:mean_fitness]=fitness_promig
    fitnessdf[:sd_fitness]=fitness_sd
    writetable(fitness_name, fitnessdf)

    effective_name=join([file_prefix, "effective.csv"],"_")
    effectivedf = DataFrame()
    effectivedf[:mean_effe]=effective
    effectivedf[:sd]=effective_sd
    writetable(effective_name, effectivedf)

    return poblacio
end

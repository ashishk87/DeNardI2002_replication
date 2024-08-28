# Loading relevant packages
using Distributions, Plots, PyPlot, LaTeXStrings, IterTools, QuantEcon, CSV, DataFrames, LinearAlgebra, Random, Statistics,Conda, CSVFiles
using Statistics, CSV, DataFrames, LinearAlgebra, Random, Statistics,Conda, BenchmarkTools, MacroTools, Base.Threads
using Optim, GLM, FreqTables, ForwardDiff, PyCall, PlotlyJS, HTTP, CovarianceMatrices, StatsBase,Printf
using JuMP, Ipopt, NLopt, StatsBase, Econometrics, CategoricalArrays, PyCall,  RCall, PrettyTables,FixedEffectModels, XLSX, BlockArrays,BlockBandedMatrices, NLsolve, Dates, IJulia, Interact
using IJulia, Base.Threads, DelimitedFiles, UnPack
using BlackBoxOptim, Interpolations
using JLD2


include("tools.jl") # loads some primitive functions that are used in the notebook




# Defining parameters for the model

# 1: Defining a parameter dictionary for the model

function params_dict(;tlength=5, dage=5, dy = 4, def = 45, TL = 9,
    dsur = 65, T = 14, TR = T-TL, 
    β = (0.95)^tlength, σ = 1.5, ϕ₁ = -9.5,ϕ₂ = 11.6, γ = 0.85, hbet = 0.67,
    τₐ = 0.2, τ_b = 0.1, 
    # exb = 8.5, 
    ggovt = 0.18, pens = 0.4,
    maxa = 50.0, mina = 0.0, da = 120, minc = 0.0, maxc = 1.0,
    σ_z = sqrt(0.3), σ_h = sqrt(0.37), 
    y = [0.2594, 0.6513, 1.5355, 3.8547],
    Qy = [0.7132 0.2764 0.0104 0.0; 0.1467 0.6268 0.2210 0.0055; 0.0055 0.2210 0.6268 0.1467; 0.0 0.0104 0.2764 0.7132],
    Qyh = [0.5257 0.4228 0.0509 0.0007; 0.1461 0.5538 0.2825 0.0176; 0.0176 0.2825 0.5538 0.1461; 0.0007 0.0509 0.4228 0.5257],
    # y = [0.3272188, 1.000000, 3.056059],
    # Qy = [0.7490493      0.2460441      4.9065556E-03;0.1650564       0.6698872  0.1650564; 4.9065561E-03    0.2460441   0.7490494],
    # Qyh = [0.6421676      0.3458525      1.1979862e-02; 0.1682615      0.6634769   0.1682615; 1.1979862E-02  0.3458525      0.6421676],
    prod_inher = true, vol_beq = true)
    """
    tlength: length of the time period
    dage: parent-child age difference
    dy: no. of income 
    def: size of eflife (which is age-efficiency profile for workers)
    TL = def/tlength: No. of working periods
    dsur = 65 # alive fro, 20 to 65
    T = 14 # Max no. of periods alive (dsur/tlength + 1)
    TR = T-TL # No. of retirement periods

    # Preferences
    σ = 1.5 # risk aversion
    β = (0.9455)^tlength # discount factor
    ϕ₁ = -9.5 # joy of bequeathing parameter 1
    ϕ₂ = 8.0 # joy of bequeathing parameter 2
    γ = income autoregressive coefficient
    hbet: inheritance autoregressive coefficient
    gpop: population growth rate 0.05773918371694153

    # Government and tech
    r: exogenous interest rate (ensures annual K/Y = 3)
    τₐ: tax rate on capital income
    τ_b: tax rate on bequests
    exb: years avg earn for estate taxes exempt
    ggovt.: government expenditure
    pensions: repl rate*avg income
    δ: depreciation rate (6% yearly)

    # Grids
    maxa: maximum level of capital
    mina: minimum level of capital
    da = dar: No. of  points on asset grid
    dc: No. of points on consumption grid
    minc: minimum cons proprtion 
    maxc: maximum cons proportion

    # Income process
    σ_z: standard deviation of income shocks
    σ_h: std deviation of inheritance shocks

    # The values for the following three objects are directly taken from the paper
    y: 4-state productivity grid (values taken from the appendix of the paper)
    Qy: transition matrix for y (values taken from the appendix of the paper)
    Qyh: inheritance transition matrices (values taken from the appendix of the paper)
    Qher: productivity transmission matrix for a 20 yr old parent to a 20 yr old child
    invher: invariant distribution of productivity at birth
    inv40: invariant distribution of productivity at 40 years of age 
    jt4020: joint distribution of parents productivity at 40 years of age and kids productivity at 20 years of age
    eff: age-efficiency profile for workers from 20 to 65 based on the elflife.dat file and tailored to then period length of 5 years
    sur: survival probabilities from age 20 to 90
    csur: cumulative survival probabilities
    numch: per capita number of children

    """
    # Set a Random seed
    Random.seed!(2022)


    gpop = (1.012)^tlength
    r = 1.06^tlength - 1
    δ = 1 - (0.94)^tlength
    dc = da


    #=
    ! matrix Qyh: inheritance of y: 40 years old parent-20 years old child
    ! want transmission  matrix Qher:
    ! 20 years old parent- her future 20 years old child
    ! (because want the invariant distribution at birth of these characterstics in population)
    ! have: c_20'=p_40'*Qyh, (i.e. children at 20 gets assigned what parents are at 40 times the inheritance transmission matrix) 
    ! We can rewrite the above as:
    ! c_20'=p_35'*Qy*Qyh, as p_40'=p_35'*Qy. (i.e. parents at 40 gets assigned what parents are at 35 times the income transmission matrix)
    ! Keep doing it recursively until we get the following
    =#
    Qher = Qy * Qy * Qy * Qy * Qyh
    # Compute the eigenvalues and eigenvectors of the transpose of Qher
    eig = eigen(Qher')
    #Extract the principal eigenvalue and eigenvector
    EVAL = eig.values
    EVEC = eig.vectors
    # Sort eigenvalues in decreasing order and rearrange eigenvectors accordingly
    sorted_indices = sortperm(EVAL, rev=true)
    EVAL = EVAL[sorted_indices]
    EVEC = EVEC[:, sorted_indices]
    # Check the conditions
    if (abs(real(EVAL[1]) - 1) > 1e-4) || (abs(real(EVAL[2] - 1)) < 1e-4)
        error("problems with Qy")
    end
    eigval = real(EVAL[1])
    eigvec = real.(EVEC[:, 1])

    # invher is the stationary distribution of productivity at age 20 when productivity inheritance is operational
    invher = eigvec ./ sum(eigvec)
    # distribution of characterstics at 40 years of age:  multiply Qy 4 times coz it gets from age 20 to 40 (20 ->25 ->30 ->35 ->40) and inverh is the invariant distribution of characterstics at age 20
    inv40 = (invher'*Qy*Qy*Qy*Qy)'


    # Calculate the joint distribution for parent's characterstics at age 40 and children's characterstics at age 20.
    #=
    So, (i, j)th entry of jt4020 is the distribution that the parent at age 40 is in the i-th income state and the child at age 20 is in the j-th income state.
    =#

    #=
    In secion 4 of the paper, the author mentions the following:

    "To be consistent across experiments, I use the same initial distribution of productivity for 20 yr old workers in all simulations.
     To do so, I compute the initial aggregate distribution of productivity implied by the experiment with inheritance of productivity 
     and use it to initialize all the 20-year-old workers in the simulations without productivity inheritance"

    invher is the initial aggregate distribution of productivity implied by the experiment with inheritance of productivity.
    - With productivity inheritance we use invher and Qy to get to the stationary distribution of productivity at age 40 (inv40 = Qy*Qy*Qy*Qy*invher).
      And then use inv40 and Qyh to get the joint distribution of parents productivity at 40 and kids productivity at 20 (jt4020).
    
    - Without productivity inheritance, we start with invher and use (Qy)⁴  to get the joint distribution of parents productivity at 40 and kids productivity at 20 (jt4020),
       rather than using inv40 and Qyh.
    =#

    jt4020 = zeros(dy, dy)
    Q_y_20_to_40 = Qy^4
    if prod_inher
        for i in 1:dy
            for j in 1:dy
                jt4020[i, j] = Qyh[i, j] * inv40[i]
            end
        end
    else
        jt4020  = inv40*invher'
    end



    #=
    Construct the conditional survival probability from age 20 to 90
    =#
    # Read survival probabilities from file
    surlife = readdlm("/Users/ashishkumar/Desktop/DeNardI2002_replication/surlife.dat")
    # Process survival probabilities
    sur = ones(T) # Initialize survival probabilities
    for i in 1:tlength:dsur
        sur[(i-1)÷tlength + 1] = prod(surlife[i:i+tlength-1])
    end

    # Adjustments for early ages
    sur[1:8] .= 1.0
    sur[end] = 0.0 

    # Cumulative survival probabilities
    csur = ones(T + 4)
    csur[5:end] = sur
    for i in 5:T+4
        csur[i] = prod(sur[1:i-4])
    end


    # age productivity transition matrix. For ex. (1, dy+1) entry of this matrix gives the prob. of moving from age
    # 1 and productivity 1 to age 2 and productivity 1, (1, dy+2) gives the prob. of moving from age 1 and productivity 1 to age 2 and productivity 2, and so on.
    EM = zeros(TL*dy,TL*dy) 
    eq, eminv, eminv1 = [zeros(TL*dy) for _ in 1:3]
    for i in 1:TL-1 # loop over today's age
        for j in 1:dy # loop over productivity
            EM[(i-1)*dy + j, dy + (i-1)*dy + 1: dy + i*dy] = sur[i]*Qy[j,:]
        end
    end
    # Now computing the invariant distribution of age-productivity
    eq[1:dy] = invher
    eminv = eq
    ϵ = 1.0
    while ϵ > 1e-8
        eminv1 = eminv
        eminv = (EM'*eminv)/gpop .+ eq
        ϵ = maximum(abs.(eminv - eminv1))
    end
    eminv = eminv./sum(eminv)

    #=
    Form a matrix with age on the rows and productivity on the columns. Every cell contains 
    the no. of people in each cell in order to compute the age-efficiency profile implied by 
    the income process.
    =#
    ageyn = zeros(TL,dy)
    ageyn .= transpose(reshape(eminv, dy, TL))
    # transform the matrix in terms of fraction of people of a given age
    for i in 1:TL 
        ageyn[i,:] = ageyn[i,:]./sum(ageyn[i,:])
    end
    # compute the avg. age productivity profile implied by income process
    endeff = zeros(TL)
    endeff = ageyn*y

    #= Construct age-efficiency profile for workers: from 20 to 65
    I use the elflife.dat file that author provided to construct this object
    =#
    # ! construct age-efficiency profile for workers
    #Reading data from a file
    # Open and read the 'eflife.dat' file
    eflife = readdlm("/Users/ashishkumar/Desktop/DeNardI2002_replication/eflife.dat")
    # Initialize the efficiency array
    eff = zeros(Float64, ceil(Int, length(eflife) / tlength))

    # Transform yearly age-efficiency profile according to the period length
    for i in 1:tlength:def
        eff[(i-1)÷tlength + 1] = sum(eflife[i:i+tlength-1])
    end

    # Normalize the efficiency profile
    eff = (eff *TL)/sum(eff)
    eff = eff ./ endeff

    #! Experiment
    # Surv prob.
    sur1, pfrac = [zeros(T)  for _ in 1:2]
    sur1[1] = 1.0
    sur1[2:T] = sur[1:T-1]
    for i in 1:T
        pfrac[i] = prod((sur1[1:i] ./ gpop))
    end

    pfrac = pfrac./sum(pfrac)
    aveinc = sum(pfrac[1:TL].*eff.*endeff)

    exb = 8*aveinc # 40 yrs of avg. earnings
    #

    # Compute per capita number of children-Same as in DeNardi (2002)
    # That is done as follows:
    # newbornsₜ₊₁ = gpop*newbornsₜ - (1) Population law of motion
    # newbornsₜ = numch*parentsₜ - (2) (numch is the per capita number of children), and therefore (2) gives us total number of children
    # parentsₜ = newborns_{t-dage}*sur[1] - (3) {No. of parents at 25 equals the number of newborns at 20 times the survival probability from 20 to 25 (one until they are 20 and sur[1] b/w 20 and 25)}
    # Using (2) and (3): newbornsₜ = numch*newborns_{t-dage}*sur[1] - (4)
    # Using (1): newbornsₜ = gpop*newbornsₜ₋₁
    #                     = (gpop)^2*newbornsₜ₋₂
    #...                  = (gpop)^(dage)*newborns_{t-dage} - (5)
    # Using (4) and (5): numch = (gpop)^dage/sur[1] - (6)
    numch = gpop^dage/sur[1]

    # Asset Grid
    agrid = range(mina, sqrt(maxa), length=da)
    agrid = agrid.^2

    # Asset Grid Net of Estate taxes
    anet = zeros(size(agrid))
    anet .= ifelse.(agrid .<= exb, agrid / numch, (agrid .- τ_b .* (agrid .- exb)) / numch)

    # Consumption Grid
    cgrid = collect(range(minc, maxc, length=dc+1))

    # Warm-Glow utility from net bequests-grid over all possible bequests levels
    ϕ = zeros(size(agrid))
    # ϕ1 = zeros(size(agrid))
    wg = zeros(size(agrid))
    if vol_beq
        # wg .= ((1 .+ (agrid .- τ_b.*max.(agrid .- exb, zeros(size(agrid))))./ ϕ₂).^(1 - σ) .- 1).*ϕ₁ 
        # ϕ1 .= [(1 + (a - τ_b * max(a - exb, 0.0)) /ϕ₂)^(1 - σ) - 1 for a in agrid] .* ϕ₁
        ϕ .=((1 .+ (agrid .- τ_b.*max.(agrid .- exb, 0.0 .* agrid))./ϕ₂).^(1-σ) .- 1).*ϕ₁
    else
        ϕ .= 0.0
    end


    return Dict(:tlength => tlength, :dage => dage, :dy => dy, :def => def, :TL => TL, :dsur => dsur, :T => T, :TR => TR,
    :β => β ,:σ => σ, :ϕ₁ => ϕ₁, :ϕ₂ => ϕ₂, :γ => γ, :hbet => hbet, :τₐ => τₐ, :τ_b => τ_b, :exb => exb, :ggovt => ggovt, :pens => pens,
    :maxa => maxa, :mina => mina, :da => da, :minc => minc, :maxc => maxc, :dc => dc, :r => r, :δ => δ, :gpop => gpop, :σ_z => σ_z, :σ_h => σ_h,
    :y => y, :Qy => Qy, :Qyh => Qyh, :Qher => Qher, :invher => invher, :inv40 => inv40, :jt4020 => jt4020, :eff => eff, :sur => sur, :csur => csur, :numch => numch,
    :agrid => agrid, :anet => anet, :cgrid => cgrid, :ϕ => ϕ)

end

params = params_dict();


# Main function that solves the model and solves the model and stores the results in a dictionary

function denardi_main(params)
    """
    This function solves and simulates the model in a general equilibrium
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params


    # Set a Random seed
    Random.seed!(2022)

    # Initializing objects such as the value function, policy functions, etc.

    #=
    distribution of bequests:
    -Bequests in this model depends on ageent's age (can inheirt upto first year of retirement)
    and on her parents productivity level at age 40.
    -Agents do not die until they are 60, so their kids can't inherit unitl they are 40.
    - beqdist contains bequest distribution gross of estate taxes
    =#
    beqdist = zeros(da, dy, TL-3) 
    beqdist[1, :, :] .= 1  # initialization

    #=
    Value function for retirees:
    -state 1: cash at hand
    -state 2: age
    =#
    vr = zeros(da, TR+1)
    # Setting the last period's value fn to the value of joy of bequeathing for that period (since agent will be dead for sure by then)
    # DeNardi (2002) assumes parents gets joy of bequeathing on net bequests
    vr[ :, TR+1] .= ϕ

    #=
    Policy function for the retirees:
    -state 1: market resources at hand
    -state 2: age
    =#
    # polr = zeros(da, TR+1)  # Check again
    polr = zeros(da, TR)

    #=
    Value function and policy fn for working people:
        -state1: cash at hand
        -state2: productivity
        -state3: parents productivity at age 40 (expectations of bequests and inheritance shape up the optimal consumption decisions)
                 (parents might be dead, so extra state added to the state space)
        -state4: age
    =#
    vt = zeros(da, dy, (dy+1), TL)
    polt = zeros(Int, da, dy, (dy+1), TL)
    # weights assigned to the left grid point next to the optimal savings choice
    w_assets = zeros(da, dy, (dy+1), TL) 
    polc = zeros(da, dy, (dy+1), TL)

    #=
    Optimal savings policy function for the retirees:
        -When agent can inherit, there's a prob. dist. on next peiod's assets linked to the relevant
        prob dist on bequests
    =#
    sopt = zeros(da, dy, (dy+1) ,TL)

    # Some other objects
    agew = zeros(da, T) # assets at age t for t = 1, 2, ..., T-1
    agewT = zeros(da) # assets at age T
    beqtax=0.0
    aveinc = 0.0
    numret = 0.0
    gbal = 0.0
    dista = zeros(da)
    avgk = 0.0
    indnoh = dy*dy*da*4 
    indret = 4*dy*dy*da + 5*dy*(dy+1)*da  
    counterbeq = 0.0


    #=
    Initializing the vector of newborns, i.e. their distribution across the state-space using the invariant distribution of productivity at birth inherited from parents.
    Since no parent die before 60, no child inherits before 40, so all newborns start off from zero wealth.

    This object changes during the experiment when there is no productivity inheritance. 
    Instead of doing this, we set qinit to be equal to the stationary distribution over income states for 20 year olds (see eminv)
    =#
    qinit = zeros(dy*dy*da)
    for i in 1:dy
        for j in 1:dy
            qinit[(i-1)*da*dy + (j-1)*da + 1] = jt4020[i, j]
        end
    end

    # Initializing invariant distribution
    #=
    # total no. of states:
      1) 4*dy*dy*da: total no. of states for working people when young and parents definitely survive next period. ∴ (#age)*(#productivity)*(#parents productivity)*(#asset)
      2) 5*dy*(dy+1)*da: total no. of states for working people when young and parents might die next period. ∴ (#age)*(#productivity)*(#parents productivity)*(#asset)
      3) 5*da: total no. of states for retirees. ∴ (#age)*(#asset)
    =#
    ninvm = 4*dy*dy*da + 5*dy*(dy+1)*da + 5*da 
    invm = zeros(ninvm)

    # Kernels to smooth the bequest distribution. Will be used later on.
    numkernel = 31
    kernelfun = zeros(da, numkernel)
    kernelreg = zeros(numkernel, da)
    kernel_smooth!(params, numkernel ,kernelfun, kernelreg)


    # Arrays to store smooth bequest distribution
    kercoeff, kercoeff2 = [zeros(numkernel, dy, TL-3) for _ in 1:2]

    #=
    Solve the value fn of retirees:
    - Retirees have already received bequests and are not going to work anymore.
    - They do not pay labor income taxes.
    - So, the state variable is just the cash at hand, and their age.
    =#
    vfi_retirees!(params, vr, polr)



    #=
    Now we fill the transition matrix (M) for the retirees:
    - M is a sparse matrix with the following structure:
    - rowindex of the transition matrix is given by rowM[counter]
    - colindex of the transition matrix is given by colM[counter]
    - value of the transition matrix is given by valM[counter]
    - counter track the states. Since optimal savings can lie between two grid points, we need to keep track of the weights assigned to the left grid point next to the optimal savings choice. So, for any given state s to any other state s', there will be two counter values attached to this one transition coz of the two grid points involved.

    here counter is is the unique identifier of each cell in the transition matrix
    =#
    
    #=
    We begin by specifying the nonzero elements in the transition matrix.
    And that is nonzeroM. Logic?

    Take the first component (da*dy*dy*3)*(dy).
    - da*dy*dy*3: total no. of states for working people when young and parents definitely survive next period. ∴ (#age)*(#productivity)*(#parents productivity)*(#asset)
    - Next period agents productivity can change, and so there are dy possible states for the next period's productivity.
    - Agents assets next period are determined by optimal policy, so it can transition to only one specific asset state.
    - Agent gets older by 1 year, so his age transitions to a deterministic next period age.
    - Since agent forms expectations of bequests based on parents productivity at age 40 which remains same even next period.
    - So, in total there are dy possible transition from any of the da*dy*dy*3 current states. 
    - Therefore, there are da*dy*dy*3*dy nonzero elements for the young working people when parents definitely survive next period.

    Similarly other components of the nonzeroM can be interpreted.
    =#
    nonzeroM =  (da*dy*dy*3)*(dy*2)+(da*dy*dy)*(dy*2+dy*da)+ 
                (da*dy*dy*4)*(dy*2+dy*da)+(da*dy*4)*(dy*2)+ 
                (da*dy*dy)*da+(da*dy)*2+(da*4)*2 

    # println(nonzeroM) # 4978080
    
    
    #=
     Starting row index for the retirees. This is the last row index of the working people.
     Can also be thought of as the total no. of states for working people.
    =#
    indret = 4*dy*dy*da + 5*dy*(dy+1)*da  

            
    rowM = zeros(Int, nonzeroM) # row index of the transition matrix
    colM = zeros(Int, nonzeroM) # column index of the transition matrix
    valM = zeros(nonzeroM) # value of the transition matrix corresponding to the row and column index

    # counter for the nonzero elements in the transition matrix.
    counter = Ref(1)

    # Update the part of transition matrix for retirees
    transition_retirees!(params, polr, indret, rowM, colM, valM, counter)


    #=
     # total no. of non zero entries in the transition matrix of the oldies
     da: the no. of states of cash-on-hand
     # And they die next period for sure, and next period savings are given by policy function, so transition to next state is deteministic.
     # Therefore, there are da nonzero elements for the oldies.
    =#
    nonzeroMT = da
    rowMT = zeros(Int, nonzeroMT) # row index of the transition matrix for oldies
    colMT = zeros(Int, nonzeroMT) # column index of the transition matrix for oldies
    valMT = zeros(nonzeroMT) # value of the transition matrix for oldies
    counter2 = Ref(1)
    # Update the transition matrix for the oldies (this is seperate matrix MT which is dfferent from M which is for the retirees and working people)
    transition_oldies!(params, polr, rowMT, colMT, valMT, counter2)

    #=
    Now we solve the value function and policy function for the working people (continuing the backward induction)
    And fill the part of the transition matrix for the working people.

    Some important points:
    - Working people can only receive bequests after they turn 40 so they form expectations of bequests and inheritance. So, we'll have to keep track of the entire bequest distribution.
    - Also the budget constraint of the govt. needs to be satisfied. So, adjust τ_l until it does.
    =#

    # τ_l = 0.19671993190545248 # initial guess for labor income tax rate
    # τ_l = 0.2823769  # initial guess for labor income tax rate
    τ_l = 0.224  # initial guess for labor income tax rate
    relax = 0.7 # relaxation parameter for updating τ_l
    ϵ₂ = 1.0 # tolerance for govt. budget constraint
    # # counter for the nonzero elements in the transition matrix.
    # counter = Ref(4*da + 1)

    while ϵ₂ > 1e-4 # loop until govt. budget constraint is satisfied
        counterbeq = 1.0

        # Iterating over bequest distribution
        ϵ₁ = 1.0 # tolerance for convergence of bequest distribution

        while (ϵ₁ > 1e-4) # loop until convergence of bequest distribution
            iter = 1
            # Store the hours, minutes and seconds before the iteration in a log file
            current_time = Time(now())
            # Extract hour, minute, and second
            IHOUR = hour(current_time)
            IMIN = minute(current_time)
            ISEC = second(current_time)

            # Open the log file in append mode
            open("time_log.txt", "a") do file
                # Write the values to the file
                write(file, "Started at $IHOUR $IMIN $ISEC\n")
            end

            "
            First we compute the value and policy function of all workers and then we will compute the transition matrix for the working people.
            We divide this task in cases.


            Case-1: Workers age TL (i.e. age 60). They are going to retire next period.
                    Currently the state space is ( parents productivity, own productivity, cash at hand).
                    Within this case there are two sub-cases.

                Case-1.1: Workers who have already inherited
                    - They have already received bequests so don't track beq dist. anymore
                
                Case-1.2: Workers who did not inherit so far but will definitely inherit next period
                    - This case requires tracking the bequest distribution

            The following function updates the value function, policy function for this case.
            "
            vfi_workers_1!(params,τ_l, vr, polt, polc, w_assets, vt, sopt, beqdist)

            "
            Now moving to next case:

            Case-2: Workers aged 55 to 35.
                - They receive income shock next period, so far that was not the case
                - They can inherit next period
                - But they cannot die next period, so they do not leave bequests next period
            "
            vfi_workers_2!(params, τ_l, polt, polc, w_assets,  vt, sopt, beqdist)

            "
            Now moving to last case:

            Case-3: Workers aged 20 to 30.
                - They receive income shocks
                - They cannot inherit next period
                - And they cannot die next period, so they do not leave bequests next period
            "
            vfi_workers_3!(params, τ_l, polt, polc, w_assets,  vt, sopt, beqdist)

            # Store the hours, minutes and seconds before the iteration in a log file
            current_time = Time(now())
            # Extract hour, minute, and second
            IHOUR = hour(current_time)
            IMIN = minute(current_time)
            ISEC = second(current_time)

            # Open the log file in append mode
            open("time_log.txt", "a") do file
                # Write the values to the file
                write(file, "val fun computed at $IHOUR $IMIN $ISEC\n")
            end

            #=
            Computed the value function and policy function for all workers. Now we need to update the transition matrix for the working people.

            Important points to note: Agents die and born each period, so rows of the transition matrix do not sum to 1 (people that die disappear from the system, and at each period a new cohort is born).
            The people's distribution over the state variables evolve accoding to:
            nₜ₊₁ = nₜ'M + qₜ₊₁ , where nₜ is the distribution of people over the state variables at time t and qₜ₊₁ is the distribution of newborns at time t+1.

            How is matrix M structured?
             - 1st variable: Agent's age
             - 2nd variable: Parent's productivity
             - 3rd variable: Agent's productivity
             - 4th variable: Agent's assets

            Before retirement, all four constitute the state space. After retirement, only age and assets constitute the state space coz they don't work and ca't inherit anymore.

            We filled in the transition matrix for the retirees and the oldies who are going to die for sure next period. Now we fill in the transition matrix for the working people now.
            =#

            rowM[4*da*2 + 1:end] .= 0
            colM[4*da*2 + 1:end] .= 0
            valM[4*da*2 + 1:end] .= 0

            # # counter for the nonzero elements in the transition matrix.
            counter = Ref(4*da*2 + 1)

            # Starting row index for the people too young to have inherited
            indnoh = dy*dy*da*4 
            
            # turn all entries in w_assets to be positive
            w_assets .= abs.(w_assets)

            "
            Case-1: Transition matrix for working people who've already inherited but are not retiring next period.
                    - In the model, the earliest an agent inherits is when he's 40 yr
            "
            transition_workers_1!(params,indnoh, polt, w_assets, rowM, colM, valM, counter)

            "
            Case-2: Transition matrix for working people that have already inherited and are retiring next period.
            "
            transition_workers_2!(params,indnoh, indret, polt, w_assets, rowM, colM, valM, counter)

            "
            Case-3: Transition matrix for working people who have not inherited so far and who do not inherit next period.
                    Age 20 to 35 included
            "
            # Case-3.1: Age 20 to 35 included
            transition_workers_3_1!(params, polt, w_assets, rowM, colM, valM, counter)
            # Case-3.2: Age 40 to 55 included
            transition_workers_3_2!(params, indnoh ,polt, w_assets, rowM, colM, valM, counter)

            "
            Case-4: Working people that have not inherited but will do so next period
            "
            # Case-4.1: Age 35 (or model age 4)
            transition_workers_4_1!(params, indnoh, sopt, rowM, colM, valM, counter, beqdist)
            # Case-4.2: Working people that inherit but are not retiring next period (age 40 to 55)
            transition_workers_4_2!(params, indnoh, sopt, rowM, colM, valM, counter, beqdist)
            # Case-4.3: Working people that inherit and are retiring next period (age 60)
            transition_workers_4_3!(params, indnoh, indret, sopt, rowM, colM, valM, counter,beqdist)


            "
            Now we have the transition matrix for the working people. We can now compute the invariant distribution of the state variables.
            "
            q = zeros(dy*dy*da)
            q .= qinit
            invm1 = similar(invm)
            invm1 .= 0.0
            err = 1.0
            while err > 1e-8
                invm1 .= invm
                invm .= 0.0

                # The following for loop is the product M'*invm. Since, M is a sparse
                # matrix, to transpose M, swe simply use colM instead of rowM.
                # invm = updated distribution
                # invm1 = old distribution
                for i in 1:nonzeroM
                    invm[colM[i]] = invm[colM[i]] + (invm1[rowM[i]]*valM[i])/gpop 
                end

                invm[1:dy*dy*da] = invm[1:dy*dy*da] +  q
                err = maximum(abs, invm - invm1)
            end

            invm .= abs.(invm) ./ sum(abs.(invm))
            q .= q./sum(q)

            # Check if invm has any NaN values
            if any(isnan, invm)
                println("invm has NaN values")
                println(sum(invm))
            end

            # Now updating the bequest distribution
            beqdist .= 0

            #=
            DeNardi(2002) assumes that children observes parents productivity at age 40 before enetering job market (i.e. at age 15).
            Therefore, the ebquest distribution will be conditional on the parents productivity at age 40: beqdist(a |yp40, t).
            To construct this object, we go through several steps:

            Step-1: Find the distribution of assets at 40 yrs of age (Since we are dealing with stationary distribution, we can use the invariant distribution of assets at 40 yrs of age.)
                    Note that the above distrution will depend on parents productivity as well as their parents productivity (kid's grandparents) at age 40.
            =#

            # The following is the joint distribution of agent's income, his parents income (also at age 40), and his assets at age 40.
            temprealv2 = zeros(da*dy*(dy+1))
            temprealv2 = invm[4*dy*dy*da + 1:4*dy*dy*da + dy*(dy+1)*da] 

            # Using the above object, we will construct the conditional distribution of assets at 40 yrs of age given parents productivity at 40 yrs of age.
            #=
            Step-1 Contd: Construct joint distribution of agents parents income at 40, and agents assets at 40 conditional on his income at 40.
                          P(yp = ypⱼ, a = aₖ|y = yᵢ) = P(yp = ypⱼ, a = aₖ, y = yᵢ)/P(y = yᵢ) - (*)
                          where:
                          P(y = yᵢ) = ∑ⱼₖ P(yp = ypⱼ, a = aₖ, y = yᵢ)
            =#
            temprealv = similar(temprealv2)
            temprealv .= 0.0

            temprealv4 = zeros(dy) # marginal (stationary) distribution of income at 40 yrs of age

            # First we construct P(y = yᵢ) by summing over all possible values of yp and a
            # To construct the marginal of agents income at 40 integrate out it from the joint distribution-temprealv2. 
            # That is what the following block is doing.
            for i in 1:dy # loop over agents income at 40
                for j in 1:dy+1 # loop over grandparents income at 40 (might be dead)
                    for k in 1:da # loop over agents assets at 40
                        temprealv4[i] = temprealv4[i] +  temprealv2[(j-1)*(dy)*da + (i-1)*da + k]
                    end
                end
            end

            # Next, we construct P(yp = ypⱼ, a = aₖ|y = yᵢ) since we already have the numerator and denominator from (*)
            for j in 1:dy+1 # loop over agents income at 40
                for k in 1:da # loop over grandparents income at 40 (might be dead)
                    for i in 1:dy # loop over agents assets at 40
                        temprealv[(j-1)*dy*da + (i-1)*da + k] = temprealv2[(j-1)*(dy)*da + (i-1)*da + k]/temprealv4[i]
                    end
                end
            end

            #=
            Step-2: Next we construct the conditional distribution of characterstics of agents at 40 conditional on their income- m(a, y40, yp40|y40)
                    Coz that is what kids use to form expectations of bequest distribution. And we iterate that object forward in time using the
                    transition matrix M, i.e. m(a', y45, yp40|y40) = m(a, y40, yp40|y40)*M
                                              m(a'', y50, yp40|y40) = m(a', y45, yp40|y40)*M and so on.
            =#
            stbeq = zeros(da*dy*(dy+1), dy) #  column i of this object will contain the conditional distribution corresponding to income level i from temprealv.   

            for i in 1:dy # loop over agents income at 40 (parents income at 40)
                for jj in 1:(dy+1) # productivity of grandparent (might be dead, so added 1 to dy)
                    stbeq[(jj-1)*da*dy .+ (i-1)*da .+ collect(1:da), i] = temprealv[(jj-1)*(dy)*da .+ (i-1)*da .+ collect(1:da)]
                end 
            end

            #=
            Now having constructed the conditional parents characterstics at age 40-m(a, y40, yp40|y40), we will use this object
            along with transition matrix to construct bequest distribution. Intuitively, we are just going to take dot product of 
            this period state with corresponding entry from transition matrix to obtain next period same distribution. And we will 
            continue to do it until we reach age 60 and beyond coz parents leave bequests after crossing age 60 (coz of mortality risk 
            they can leave bequest anytime after that).
            =#
            beq_dist_update_1!(params, indnoh, rowM, colM, valM, stbeq)

            #=
            So, now we have updated the conditional distribution of parents characterstics m(a, y40, yp40|y40) upto age 60.
            That's what previous block of the code was all about. Now we iterate this distribution one step ahead, i.e. when agents attain age 65.
            So, like the above block we consider two cases:
              (i) working that have already inherited
              (ii) people that have not yet inherited (but will inherit next period since everyone receive bequests by the time they retire)
            =#
            stbeqret = zeros(da, dy)
            beq_dist_update_2!(params, indnoh, indret, nonzeroM, rowM, colM, valM, stbeq, stbeqret)

            # The updated stbeqret serves as the bequest distribution if agents die aged 65 (right when they retire)
            beqdist[:, :, TL-dage-3] .= stbeqret

            #=
            So far, we've computed asset distribution of agents until they turn 65 (when they start facing mortality risk).
             Now we will further rollover this asset distribution further (age 70, 75, 80, 85) using relevant part of transition matrix.
            =#
            beq_dist_update_3!(params, indret, nonzeroMT, rowM, colM, valM, rowMT, colMT, valMT, stbeqret, beqdist)

            # Smooth bequest distribution with Kernels
            kercoeff2 .= kercoeff
            beqdist2 = copy(beqdist)
            kercoeff .= 0.0
            for i in 1:TL-3
                kercoeff[ :, :, i] = kernelreg*beqdist[ :, :, i]
                beqdist[ :, :, i] = kernelfun*kercoeff[ :, :, i] 
            end

            # Set negative entries to zero
            beqdist[beqdist .< 0] .= 0.0

            # Normalize smoothed bequest distribution
            for i in 1:TL-3
                for jj in 1:dy
                    beqdist[:, jj, i] .= beqdist[:, jj, i]./sum(beqdist[:, jj, i])
                end
            end

            # # Save kercoeff and kercoeff2 in matrices
            # @save "kercoeff.jld2" kercoeff
            # @save "kercoeff2.jld2" kercoeff2

            ϵ₁ = maximum(abs, kercoeff - kercoeff2)
            @info "ϵ₁: $ϵ₁"
            @info "counterbeq: $counterbeq"
            counterbeq += 1
        end # end of while loop for convergence of bequest distribution

        # Store the hours, minutes and seconds before the iteration in a log file
        current_time = Time(now())
        # Extract hour, minute, and second
        IHOUR = hour(current_time)
        IMIN = minute(current_time)
        ISEC = second(current_time)

        # Open the log file in append mode
        open("time_log.txt", "a") do file
            # Write the values to the file
            write(file, "beqdist computed at $IHOUR $IMIN $ISEC\n")
        end

        #=
        Compute the govt. budget constraint
        =#

        # First compute distribution over assets
        dista .= 0.0
        for i in 1:4*dy*dy + 5*(dy+1)*dy + 5  # sum people up
            dista .= dista .+ invm[(i-1)*da + 1:i*da]
        end

        # Avg k: multiply each asset level for no. of people owning it, working or retired
        avgk = sum(dista .* agrid)

        #=
        ! construct a matrix with assets on the rows and age on the columns
        !        age 1  age 2 ......age TL.....ageT
        !  a1
        !  a2
        !  ...
        !  ada
        =#
        agew = zeros(da, T) # assets at age t for t = 1, 2, ..., T-1
        agewT = zeros(da) # assets at age T
        indnoh = dy*dy*da*4
        indret = 4*dy*dy*da + 5*dy*(dy+1)*da
        asset_by_age!(params, invm, indnoh, indret, nonzeroMT, rowMT, colMT, valMT, agew, agewT)

        # Count no. of assets in agrid less than exb
        tempint1 = sum(agrid .< exb)
        beqtax = 0.0
        for i in tempint1+1:da 
            tempreal = τ_b*(agrid[i] - exb)
            for j in TL+1:T
                beqtax += tempreal * agew[i, j] * (1 - sur[j-1]) / sur[j-1]
            end
            # 90 yrs old
            beqtax += tempreal * agewT[i]
        end

        # Avg. income and num retirees
        aveinc, numret, pfrac = age_prod_inv_dist(params)

        # Now calculating govt. balance
        gbal = τₐ*r*avgk + τ_l*aveinc - pens*numret - ggovt*(aveinc + r*avgk) + beqtax
        # @info "avek: $avgk"
        # @info "gbal: $gbal"
        τ_l_2 = τ_l
        τ_l = τ_l - relax*(gbal/aveinc)
        ϵ₂ = abs(τ_l - τ_l_2)
        @info "ϵ₂: $ϵ₂"
        @info "τ_l: $τ_l"
        @info "τ_l_2: $τ_l_2"
        @info "gbal: $gbal"
        @info "aveinc: $aveinc"


        
         
    end # end of while loop for govt. budget constraint

    # Store the hours, minutes and seconds before the iteration in a log file
    current_time = Time(now())
    # Extract hour, minute, and second
    IHOUR = hour(current_time)
    IMIN = minute(current_time)
    ISEC = second(current_time)

    # Open the log file in append mode
    open("time_log.txt", "a") do file
        # Write the values to the file
        write(file, "New  Tau_l computed at $IHOUR $IMIN $ISEC\n")
    end 

    # Return every object that can be used for further analysis
    return Dict(:vr => vr, :polr => polr, :polt => polt, :polc => polc, :w_assets => w_assets, :vt => vt, :sopt => sopt, :beqdist => beqdist,  :invm => invm, :ninvm => ninvm, :kercoeff => kercoeff,  :agew => agew, :agewT => agewT, :dista => dista, :avgk => avgk, :aveinc => aveinc, :numret => numret, :gbal => gbal, :τ_l => τ_l, :beqtax => beqtax,
    :indret => indret, :indnoh => indnoh, :rowM => rowM, :colM => colM, :valM => valM, :counter => counter, :rowMT => rowMT, :colMT => colMT, :valMT => valMT, :counter2 => counter2, :counterbeq => counterbeq)

end

"
First preparing results for the model with both links: Parents bequest motive and productivity inheritance.
"

# Run the main function and store the results. The output of the function is a dictionary containing all the objects that can be used for further analysis.

#1) Both intergenrational links: Parents bequest motive and productivity inheritance

params = params_dict()
# @save "params_dict.jld2" params

#2251.081121 seconds (61.33 G allocations: 1.500 TiB, 3.48% gc time, 0.20% compilation time: 1% of which was recompilation)
#5018.040279 seconds (148.58 G allocations: 3.634 TiB, 4.77% gc time, 0.01% compilation time)
@time results = denardi_main(params)
# @save "baseline_model_dict_matched.jld2" results
# @load "baseline_model_dict_matched.jld2" results
wgini, wealth_share = wealth_gini(params, results; exclude_young = true)
twr, trwealth, avgk = transfer_wealth_ratio(params, results)


#2)One link:Parents bequest motive but no productivity inheritance
params_exp_1 = params_dict(;prod_inher = false, β = (0.95)^5)
results_exp_1 = denardi_main(params_exp_1)
# @save "No_prod_inher_model_dict.jld2" results_exp_1
# @load "No_prod_inher_model_dict.jld2" results_exp_1
wgini1, wealth_share1 =   wealth_gini(params_exp_1, results_exp_1; exclude_young = true)
twr1, trwealth1, avgk1 = transfer_wealth_ratio(params_exp_1, results_exp_1)

#3) One link: Productivity inheritance and no bequest motive
params_exp_2 = params_dict(;β = (0.96)^5, vol_beq = false)
results_exp_2 = denardi_main(params_exp_2)
# @save "No_beq_motive_model_dict.jld2" results_exp_2
# @load "No_beq_motive_model_dict.jld2" results_exp_2
wgini2, wealth_share2 =   wealth_gini(params_exp_2, results_exp_2; exclude_young = true)
twr2, trwealth2, avgk2 = transfer_wealth_ratio(params_exp_2, results_exp_2)

#4) No intergenrational links: No bequest motive and no productivity inheritance-unequal (accidental) bequests to children
params_exp_3 = params_dict(;β = (0.96)^5,  vol_beq = false, prod_inher = false)
results_exp_3 = denardi_main(params_exp_3)
# @save "No_beq_motive_no_prod_inher_model_dict.jld2" results_exp_3
# @load "No_beq_motive_no_prod_inher_model_dict.jld2" results_exp_3
wgini3, wealth_share3 =   wealth_gini(params_exp_3, results_exp_3; exclude_young = true)
twr3, trwealth3, avgk3 = transfer_wealth_ratio(params_exp_3, results_exp_3)
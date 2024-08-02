# Loading relevant packages
using Distributions, Plots, PyPlot, LaTeXStrings, IterTools, QuantEcon, CSV, DataFrames, LinearAlgebra, Random, Statistics,Conda, CSVFiles
using Statistics, CSV, DataFrames, LinearAlgebra, Random, Statistics,Conda, BenchmarkTools, MacroTools, Base.Threads
using Optim, GLM, FreqTables, ForwardDiff, PyCall, PlotlyJS, HTTP, CovarianceMatrices, StatsBase,Printf
using JuMP, Ipopt, NLopt, StatsBase, Econometrics, CategoricalArrays, PyCall,  RCall, PrettyTables,FixedEffectModels, XLSX, BlockArrays,BlockBandedMatrices, NLsolve, Dates, IJulia, Interact
using IJulia, Base.Threads, DelimitedFiles, UnPack
using BlackBoxOptim, Interpolations





# Define the function binary_search with types specified for each input
function binary_search(x::Vector{Float64} , xi::Float64)
    """
    inputs:
    x: Grid on which search is to be performed
    xi: Point for which we want to find the closest grid point

    output:
    ix: Index of the closest grid point to the left of xi

    """
    Nx = length(x)
    
    # a. checks
    # Check if the target value (xi) is less than or equal to the first element in the array
    # If so, return 0 as it's the first position
    if xi <= x[1]
        return 1
    # Check if the target value (xi) is greater than or equal to the second last element in the array
    # If so, return Nx-2 as it's the position before the last
    elseif xi >= x[Nx-1]
        return Nx-1
    end
    
    # b. binary search
    # Initialize half to half of Nx, this will be used to reduce the search space
    half = Nx ÷ 2  # Note: ÷ is integer division in Julia, equivalent to // in Python
    
    imin = 1  # Initialize the lower bound to 1
    
    while half > 0
        imid = imin + half  # Calculate the middle index
        # If the middle element is less than or equal to the target, move the search to the right half
        if x[imid] <= xi
            imin = imid
        end
        # Reduce the search space by half
        Nx -= half
        half = Nx ÷ 2
    end
    
    return imin
end

function kernel_smooth!(params, numkernel ,kernelfun, kernelreg)
    """
     Implementing a kernel smoothing technique to smooth out the 
     distribution of bequests (beqdist) over an asset grid.kernelfun and kernelreg will 
     be used later on for that purpose.

     Output: (objects to be modified)
        - kernelfun
        - kernelreg
    """

    # Unpacking parameters 
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    # Parameters for the kernel smoothing
    # numkernel = 31 # number of kernel functions
    groupoint = da/(numkernel - 1) # no. of points b/w the mean of 1 kernel and mean of next
    middlepoint = floor(Int, (groupoint + 1)/2) # middle point of the kernel
    numstdev = 0.5 

    indmeank, upint, lowint = zeros(numkernel-1), zeros(numkernel-1), zeros(numkernel-1)

    # Step-1: Find location of mean and each upper and lower bound of reference intervals
    for i in 1:numkernel-1 
        indmeank[i] = middlepoint + (i-1)*groupoint
        upint[i] = groupoint + (i-1)*groupoint
        lowint[i] = 1 + (i-1)*groupoint
    end

    meankernel, stdevk = [zeros(numkernel-1) for _ in 1:2]
    meankernel .= agrid[indmeank]
    stdevk .= (agrid[upint] - agrid[lowint])/2*numstdev
    kernelfun .= 0.0
    kernelfun[1,1] = 1.0 # first kernel is degenerate on zero

    # Approximate the normal dist. that constitute the kernel on our grid
    for j in 1:numkernel-1
        
        kernelfun[i, j+1] =  ((agrid[2] - meankernel[j]) / (agrid[2] - agrid[1])) * 
                             (cdf(Normal(0, 1), (agrid[2] - meankernel[j]) / stdevk[j]) - 
                             cdf(Normal(0, 1), (a[1] - meankernel[j]) / stdevk[j])) + 
                            (stdevk[j] / (sqrt(2.0 * pi) * (agrid[2] - agrid[1]))) * 
                            (exp(-((agrid[2] - meankernel[j])^2) / (2.0 * stdevk[j]^2)) - 
                             exp(-((agrid[1] - meankernel[j])^2) / (2.0 * stdevk[j]^2)))

        for i in 2:da-1
            kernelfun[i, j+1] = ((meankernel[j] - agrid[i-1]) / (agrid[i] - agrid[i-1])) * 
                                (cdf(Normal(0, 1), (agrid[i] - meankernel[j]) / stdevk[j]) - 
                                cdf(Normal(0, 1), (agrid[i-1] - meankernel[j]) / stdevk[j])) - 
                                (stdevk[j] / (sqrt(2 * π) * (agrid[i] - agrid[i-1]))) * 
                                (exp(-((agrid[i] - meankernel[j])^2) / (2 * stdevk[j]^2)) - 
                                exp(-((agrid[i-1] - meankernel[j])^2) / (2 * stdevk[j]^2))) + 
                                ((agrid[i+1] - meankernel[j]) / (agrid[i+1] - agrid[i])) * 
                                (cdf(Normal(0, 1), (agrid[i+1] - meankernel[j]) / stdevk[j]) - 
                                cdf(Normal(0, 1), (agrid[i] - meankernel[j]) / stdevk[j])) + 
                                (stdevk[j] / (sqrt(2 * π) * (agrid[i+1] - agrid[i]))) * 
                                (exp(-((agrid[i+1] - meankernel[j])^2) / (2 * stdevk[j]^2)) - 
                                exp(-((agrid[i] - meankernel[j])^2) / (2 * stdevk[j]^2)))

        end

        kernelfun[da, j+1] = ((meankernel[j] - agrid[da-1]) / (agrid[da] - agrid[da-1])) *
                             (cdf(Normal(0, 1), (agrid[da] - meankernel[j]) / stdevk[j]) - 
                             cdf(Normal(0, 1), (agrid[da-1] - meankernel[j]) / stdevk[j])) - 
                             (stdevk[j] / (sqrt(2.0 * π) * (agrid[da] - agrid[da-1]))) * 
                             (exp(-((agrid[da] - meankernel[j])^2) / (2.0 * stdevk[j]^2)) - 
                             exp(-((agrid[da-1] - meankernel[j])^2) / (2.0 * stdevk[j]^2)))

    end

    kertemp = transpose(kernelfun) * kernelfun     
    kertemp_inv = inv(kertemp)                       
    kernelreg .= kertemp_inv * transpose(kernelfun)    

end


function age_prod_inv_dist(params)
    """
    Function to compute the age-productivity inverse distribution

    Inputs:
    - params: dictionary of parameters

    Output:
    - aveinc: average income
    - numret: no. of retirees
    """

    # Unpacking parameters 
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    # age productivity transition matrix. For ex. (1, dy+1) entry of this matrix gives the prob. of moving from age
    # 1 and productivity 1 to age 2 and productivity 1, (1, dy+2) gives the prob. of moving from age 1 and productivity 1 to age 2 and productivity 2, and so on.
    EM = zeros(TL*dy,TL*dy) 
    eq, eminv, eminv1 = [zeros(TL*dy) for _ in 1:3]
    for i in 1:TL-1 # loop over today's age
        for j in 1:dy # loop over productivity
            EM[(i-1)*dy + j, dy + (i-1)*dy + 1: dy + i:dy] = sur[i]*Qy[j,:]
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

    # Normalizing eff
    eff = eff./endeff

    # Surv prob.
    sur1, pfrac = [zeros(T)  for _ in 1:2]
    sur1[1] = 1.0
    sur[2:T] = sur[1:T-1]
    for i in 1:T
        pfrac[i] = prod(sur1[1:i] ./ gpop)
    end
    pfrac = pfrac./sum(pfrac)
    aveinc = sum(pfrac[1:TL].*eff.*endeff)
    numret = sum(pfrac[TL+1:T])

    return aveinc, numret
end

u(c,σ) = (c^(1-σ))/(1-σ)


function interplin(x, y, z)
    # Create a linear interpolation object
    itp = interpolate((x,), y, Gridded(Linear()))
    
    # Create an extrapolation object with linear extrapolation
    eitp = extrapolate(itp, Line())

    # Use the extrapolation object to get interpolated/extrapolated values
    v = [eitp(zi) for zi in z]
    
    return v
end


function vfi_retirees!(params, vr, polr)
    """
    Inputs:
    - params: dictionary of parameters
    - vr: value function of retirees (which will be updated)
    - polr: policy function of retirees (which will be updated)

    Outputs:
    - v: value function of retirees

    Note: -Retirees do not work and don't pay labr income taxes
    """
    # Unpacking parameters 
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params
    
    #=
    Different choices of consumption tmrw
    =#
    cr = zeros(da, da)
    for i in 1:da # looping over cash-on hand states
        cr[ :, i] = (1 + r*(1 - τₐ))*agrid[i] .+ pens .- agrid # cr[j, i] here gives the consumption if past savings are agrid[i], and this period's savings are agrid[j]
    end

    # taking care of negative consumption and defining utilities
    ur = zeros(da, da)
    ur[cr .> 0] .= u.(cr[cr .> 0], σ)
    ur[cr .<= 0] .= -1e10

    valr = zeros(da)

    # Backward induction step
    for j in TR:-1:1 # loop over age
        for i in 1:da   # loop over cash-on-hand
            valr .= ur[:,i] .+ sur[TL+j]*β*vr[ :, j+1] + (1 - sur[TL+j]).*ϕ

            # Choose the maximum index of the value function
            imax = argmax(valr)
            polr[i,j] = imax  # Here polr gives the optimal index of the asse choice coz cr[ :,i] and hence ur[:,i] are defined over the agrid, so the maximum element says if the savings today is ar[j] then cr[j, i] is the optimal consumption
            vr[i,j] = valr[imax]     # value function: maximum value
        end
    end
end

function transition_retirees!(params, polr, indret, rowM, colM, valM, counter)
    """
    Function to compute the part of transition matrix of retirees

    Inputs:
    - params: dictionary of parameters
    - indexr: policy function of retirees. Gives the index of the asset choice of retirees as a function of their cash-on-hand and age.
    - indret: tital no. of states when agent was young and working
    

    Output- This function modifies the following objects:
    - rowM: add row indices of the nonzero elements of the transition matrix of retirees
    - colM: add column indices of the nonzero elements of the transition matrix of retirees
    - valM: add values of the nonzero elements of the transition matrix of retirees
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of retirees
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    # Step-1: Compute the transition matrix of retirees

    for j in 1:TR-1 # loop over age
        for i in 1:da # loop over cash-on-hand

            #=
             Logic for rowM:
                -indret is the total no. of states when agent was young and working. So we need to add this coz we are in the part of retirees' transition matrix
                - i is index of the cash-on-hand state of retirees
                - For any given age j, there are da states of cash-on-hand. So, when we increment age by 1, we need to add da to the row index
            =#
            rowM[counter] = indret + (j - 1)*da + i
            #=
             Logic for colM:
                -Same logic as above. The only difference is:
                   1) tmrw's asset choice comes from policy function, so we use the indexr matrix
                   2) Since agent's 1 year older now, we ignore the first da states of cash-on-hand corresponding to the agents forgone age. Those elements of the transition matrix are always 0.
            =#
            colM[counter] = indret +  j*da  + polr[i,j]
            #=
            Logic for valM:
                - To get to the next state from their existing state, the only uncertainity they face concerns their survival probability.
            =#
            valM[counter] = sur[j + TL]
            counter += 1
        end

    end
    
end


function transition_oldies!(params, polr, rowMT, colMT, valMT, counter2)
    """
    Function to compute the part of transition matrix of oldies (who will die next period for sure)

    Inputs:
    - params: dictionary of parameters
    - polr: policy function of retirees. Gives the index of the asset choice of retirees as a function of their cash-on-hand and age.
    

    Output- This function modifies the following objects:
    - rowMT: add row indices of the nonzero elements of the transition matrix of oldies
    - colMT: add column indices of the nonzero elements of the transition matrix of oldies
    - valMT: add values of the nonzero elements of the transition matrix of oldies
    - counter2: counter to keep track of the number of nonzero elements of the transition matrix of oldies
    """
    for jjj in 1:da # asset state of oldies
        rowMT[counter2] = jjj # row index 
        colMT[counter2] = polr[jjj,T-TL] # column index
        valMT[counter2] = 1  # value
        counter2 += 1
    end

end


function vfi_workers_1!(params,τ_l, vr, polt, polc, w_assets, vt, sopt, beqdist)
    """
    Function to compute the value and policy function of workers aged TL
    
    Inputs:
    - params: dictionary of parameters
    - τ_l: guess of labor income tax rate 
    - vr: value function of retirees (reqd for the step of backward induction)
    - beqdist: bequest distribution matrix

    Outputs- The following objects are updated:
    - vt: value function of workers aged TL-1
    - polt: policy function of workers (optimal asset index) aged TL-1
    - polc: policy function of workers (optimal consumption index) aged TL-1
    - sopt: optimal asset choice of workers aged TL-1
    - w_assets: weights assigned to two asset grid points on agrid surrounding the optimal savings
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    
    for jjj in 1:dy # loop over agents current prod. state

        # Calculate net earnings
        netearn = (1 - τ_l)*eff[TL]*y[jjj]

        for j2 in 1:da # loop over cash-on-hand of agent

            "
            Case-1: Agent has already inherited in the past
            "
            # Compute today's wealth (or market resources)
            weal = (1 + (1 - τₐ)*r)*agrid[j2] + netearn
            # Compute utility matrix for each possinle level of consumption
            c, ut, st, at, vst, phit = (zeros(dc) for _ in 1:6) # Actual consumption values as cgrid is the proportion of wealth consumed
            valr = zeros(da)
            c = cgrid[2:end]*weal
            ut = u.(c, σ)
            # Now when agent is not inheriting (inherited already), compute tmrw's assets (=savinngs) for each possible level of consumption today
            st = weal - c
            # The above savings might not fall onto the defined agrid, so we need to interpolate/extrapolate the value of savings based on the value function of retirees (backward induction step)
            vst = interplin(agrid, vr[:,1], st)
            # In case agent dies, then their savings are inherited by their children. Calculate the warm-glow utility they receive from that (net of taxes)
            phit .= (  (1 .+ (st .- τ_b.*max.(st .- exb, zeros(size(st)))) ./ ϕ₂).^(1 - σ)   .- 1).*ϕ₁ 
            # Value function
            valr .= ut .+ β.*sur[TL].*vst .+ (1 - sur[TL]).*phit
            imax = argmax(valr) # Choose the maximum index of the value function
            polc[j2, jjj, dy+1, TL] = imax
            vt[j2, jjj, dy+1, TL] = valr[imax] # value function: maximum value
            sopt[j2, jjj, dy+1, TL] = (1 - (imax/dc))*weal # optimal asset choice
            polt[j2, jjj, dy+1, TL] = binary_search(agrid, Float64(sopt[j2, jjj, dy+1, TL])) # optimal asset choice index

            "
            Case-2: Agent has not inherited in the past and will surely inherit next period
            "
            # We compute tmrw's assets for each possible level of bequests received tmrw (same grid as assets) 
            # and each possible level of consumption today
            for jj in 1:dy
                evat = zeros(dc)
                for j1 in 1:da # each possible bequest level, net of tax 
                    at .= st .+ anet[j1]./csur[TL+4]
                    vat = interplin(agrid, vr[:,1], at)
                    # Agents inherit next period, which is the last period they could have inherited, meaning their parents are at the end of their life next period. 
                    evat = evat .+ beqdist[j1, jj, TL-3].*vat # expected value of assets tmrw
                end # loop over bequests end

                valr .= ut .+ β.*sur[TL].*evat .+ (1 - sur[TL]).*phit
                imax = argmax(valr) # Choose the maximum index of the value function
                polc[j2, jjj, jj, TL] = imax
                vt[j2, jjj, jj, TL] = valr[imax] # value function: maximum value
                sopt[j2, jjj, jj, TL] = (1 - (imax/dc))*weal # optimal asset choice
                polt[j2, jjj, jj, TL] = binary_search(agrid, Float64(sopt[j2, jjj, jj, TL])) # optimal asset choice index
                w_assets[j2, jjj, jj, TL] = (sopt[j2, jjj, jj, TL] - agrid[polt[j2, jjj, jj, TL]])/(agrid[polt[j2, jjj, jj, TL]+1] - agrid[polt[j2, jjj, jj, TL]]) # weight assigned to the two asset grid points on agrid surrounding the optimal savings
            end # loop over parents income ends

        end # loop over cash-on-hand ends

    end # loop over own productivity ends

end # function ends


function vfi_workers_2!(params, τ_l, polt, polc,w_assets,  vt, sopt, beqdist)
    """
    Function to compute the value and policy function of workers aged TL-1 (55) to 4 (35)
    
    Inputs:
    - params: dictionary of parameters
    - τ_l: guess of labor income tax rate 
    - vr: value function of retirees (reqd for the step of backward induction)
    - beqdist: bequest distribution matrix

    Outputs- The following objects are updated:
    - vt: value function of workers aged TL-1
    - polt: policy function of workers (optimal asset index) aged TL-1
    - polc: policy function of workers (optimal consumption index) aged TL-1
    - sopt: optimal asset choice of workers aged TL-1
    - w_assets: weights assigned to two asset grid points on agrid surrounding the optimal savings
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

  
    for j0 in TL-1:-1:4 # loop over age
        for jjj in 1:dy # loop over agents current prod. state
            "
            Case-2.1: Agent has already inherited in the past
                      - So, parents productivity no longer matters
                      - But their own productivity next period matters for tmrws value (to compute the expected value tmrw)
            "
            vti = zeros(da)
            # vt comes from the backward induction step of workers aged TL (see last fn)
            # vt for other ages will be updated below
            vti = vt[ :, :, dy+1, j0+1]*Qy[jjj, :] # Expected value tmrw from today's perspective if current productivity is y[jjj]

            "
            Case-2.2: Agent has not inherited in the past
                      - So, they might or might not inherit next period
                      - We account for both the possibilities
                      - So, parents productivity matters
                      - And their own productivity next period matters for tmrws value (to compute the expected value tmrw)
            "
            vtai = zeros(da, dy)
            for jj in 1:dy # loop over parents productivity
                vtai[:, jj] = vt[ :, :, jj, j0+1]*Qy[jj, :] # Expected value tmrw from today's perspective if parents productivity is y[jj]
            end

            netearn = (1 - τ_l)*eff[j0]*y[jjj]

            for j2 in 1:da # looping over cash-on-hand state variable
                weal = (1 + (1 - τₐ)*r)*agrid[j2] + netearn

                # Compute utility matrix for each possinle level of consumption
                c, ut, st, at, vst, vat, valr = (zeros(dc) for _ in 1:7) # Actual consumption values as cgrid is the proportion of wealth consumed
                c = cgrid[2:end]*weal
                ut = u.(c, σ)

                "
                For case- 2.1, compute tmrw's savings for each possible level of consumption today
                "
                st .= weal .- c
                # Interpolate tmmrw's value fn over the implied savings using the value fn over the agrid
                vst = interplin(agrid, vti, st)
                # Value function today
                valr .= ut .+ β.*vst
                imax = argmax(valr) # Choose the maximum index of the value function
                polc[j2, jjj, dy+1, j0] = imax
                vt[j2, jjj, dy+1, j0] = valr[imax] # value function: maximum value
                sopt[j2, jjj, dy+1, j0] = (1 - (imax/dc))*weal # optimal asset choice
                polt[j2, jjj, dy+1, j0] = binary_search(agrid, Float64(sopt[j2, jjj, dy+1, j0])) # optimal asset choice index
                w_assets[j2, jjj, dy+1, j0] = (sopt[j2, jjj, dy+1, j0] - agrid[polt[j2, jjj, dy+1, j0]])/(agrid[polt[j2, jjj, dy+1, j0]+1] - agrid[polt[j2, jjj, dy+1, j0]]) # weight assigned to the two asset grid points on agrid surrounding the optimal savings

                "
                For Case-2.2, we need to incorporate the expected value of assets tmrw based on the parents productivity
                "
                for jj in 1:dy # loop over parents productivity
                    evat = zeros(dc)
                    for j1 in 1:da # each possible bequest level, net of tax 
                        at .= st .+ anet[j1]./csur[j0+4]
                        vat = interplin(agrid, vti[:, jj], at)
                        evat = evat .+ beqdist[j1, jj, j0-3].*vat # expected value if bequests are received tmrw
                    end # loop over bequests end
                    vst .= 0.0
                    valr .= 0.0
                    vst  = interplin(agrid, vtai[:, jj], st) # Expected value if bequests are not received tmrw
                    valr .= ut .+ β.sur[j0 + dage].*vst .+ (1 - sur[j0 + dage]).*β.*evat
                    imax = argmax(valr) # Choose the maximum index of the value function
                    polc[j2, jjj, jj, j0] = imax
                    vt[j2, jjj, jj, j0] = valr[imax] # value function: maximum value
                    sopt[j2, jjj, jj, j0] = (1 - (imax/dc))*weal # optimal asset choice
                    polt[j2, jjj, jj, j0] = binary_search(agrid, Float64(sopt[j2, jjj, jj, j0])) # optimal asset choice index
                    w_assets[j2, jjj, jj, j0] = (sopt[j2, jjj, jj, j0] - agrid[polt[j2, jjj, jj, j0]])/(agrid[polt[j2, jjj, jj, j0]+1] - agrid[polt[j2, jjj, jj, j0]]) # weight assigned to the two asset grid points on agrid surrounding the optimal savings
                end # loop over parents productivity ends
            end # loop over cash-on-hand ends
        end # loop over own productivity ends
    end # loop over age ends
end # function ends

function vfi_workers_3!(params, τ_l, polt, polc, w_assets,  vt, sopt, beqdist)
    """
    Function to compute the value and policy function of workers aged 3 (30) to 1 (20)
    
    Inputs:
    - params: dictionary of parameters
    - τ_l: guess of labor income tax rate 
    - vr: value function of retirees (reqd for the step of backward induction)
    - beqdist: bequest distribution matrix

    Outputs- The following objects are updated:
    - vt: value function of workers 
    - polt: policy function of workers (optimal asset index)
    - polc: policy function of workers (optimal consumption index)
    - w_assets: weights assigned to two asset grid points on agrid surrounding the optimal savings
    - sopt: optimal asset choice of workers 
    """

    # Unpacking parameters
       
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    for j0 in 3:-1:1 # loop over age
        for jjj in 1:dy # loop over agent's productivity
            netearn = (1 - τ_l)*eff[j0]*y[jjj]
            for j2 in 1:da # loop over cash-on-hand
                weal = (1 + (1 - τₐ)*r)*agrid[j2] + netearn

                # Compute utility matrix for each possinle level of consumption
                c, ut, st, vti, vst, valr = (zeros(dc) for _ in 1:6) # Actual consumption values as cgrid is the proportion of wealth consumed
                c = cgrid[2:end]*weal
                ut = u.(c, σ)
                # Compute tmrw's savings for each possible level of consumption today
                st .= weal .- c
                for jj in 1:dy # parents productivity
                    vti .= 0.0
                    vti = vt[ :, :, jj, j0+1]*Qy[jjj,: ] # Expected value tmrw from today's perspective if current productivity is y[jjj]
                    # Interpolate value fn (do not inherit)
                    vst = interplin(agrid, vti, st)
                    # Young value function
                    valr .= ut .+ β.*vst
                    imax = argmax(valr) # Choose the maximum index of the value function
                    polc[j2, jjj, jj, j0] = imax
                    vt[j2, jjj, jj, j0] = valr[imax] # value function: maximum value
                    sopt[j2, jjj, jj, j0] = (1 - (imax/dc))*weal # optimal asset choice
                    polt[j2, jjj, jj, j0] = binary_search(agrid, Float64(sopt[j2, jjj, jj, j0])) # optimal asset choice index
                    w_assets[j2, jjj, jj, j0] = (sopt[j2, jjj, jj, j0] - agrid[polt[j2, jjj, jj, j0]])/(agrid[polt[j2, jjj, jj, j0]+1] - agrid[polt[j2, jjj, jj, j0]]) # weight assigned to the two asset grid points on agrid surrounding the optimal savings
                end # loop over parents productivity ends
            end # loop over cash-on-hand ends
        end # loop over productivity ends
    end # loop over age ends
end # function ends


function transition_workers_1!(params,indnoh, polt, w_assets, rowM, colM, valM, counter)
    """
    Function to fill in the transition matrix for the young who've already inherited but not retiring next period
    (i.e. aged 40 to 55)
    
    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - polt: policy function of workers (optimal asset index)
    - w_assets: weights assigned to two asset grid points on agrid surrounding the optimal savings

    Outputs- The following objects are updated:
    - rowM: add row indices of the nonzero elements of the transition matrix of workers
    - colM: add column indices of the nonzero elements of the transition matrix of workers
    - valM: add values of the nonzero elements of the transition matrix of workers
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of workers
    
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    polt[polt .>= da] .= da - 1


    for jj in 5:TL-1 # loop over age
        for jjj in 1:dy # loop over agent's productivity
            for j4 in 1:da # loop over cash-on-hand
                for j5 in 1:dy # loop over agent's productivity tmrw
                    rowM[counter] = indnoh + (jj - 5)*dy*(dy+1)*da + dy*dy*da + (jjj - 1)*dy*da + j4
                    colM[counter] = indnoh + (jj - 4)*dy*(dy+1)*da + dy*dy*da + polt[j4, jjj, dy+1, jj]
                    valM[counter] = Qy[jjj, j5]*w_assets[j4, jjj, dy+1 , jj]
                    counter += 1

                    rowM[counter] = indnoh + (jj - 5)*dy*(dy+1)*da + dy*dy*da + (jjj - 1)*dy*da + j4
                    colM[counter] = indnoh + (jj - 4)*dy*(dy+1)*da + (j5 - 1)*dy*da + (polt[j4, jjj, dy+1, jj] + 1)
                    valM[counter] = Qy[jjj, j5]*(1 - w_assets[j4, jjj, dy+1 , jj])
                    counter += 1
                end
            end
        end 
        
    end
end


function transition_workers_2!(params,indnoh, indret, polt, w_assets, rowM, colM, valM, counter)
    """
    Transition matrix for working people that have already inherited and are retiring next period.
    (i.e. aged 60)
    
    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - indret: total no. of states when agent was young and working
    - polt: policy function of workers (optimal asset index)
    - w_assets: weights assigned to two asset grid points on agrid surrounding the optimal savings

    Outputs- The following objects are updated:
    - rowM: add row indices of the nonzero elements of the transition matrix of workers
    - colM: add column indices of the nonzero elements of the transition matrix of workers
    - valM: add values of the nonzero elements of the transition matrix of workers
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of workers
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    polt[polt .>= da] .= da - 1

    for jjj in 1:dy # loop over agents current prod. state
        for j4 in 1:da # loop over cash-on-hand

            rowM[counter] = indnoh + (TL-5)*dy*(dy+1)*da + dy*dy*da + (jjj - 1)*da + j4
            colM[counter] = indnoh +  polt[j4, jjj, dy+1, TL]
            valM[counter] = sur[TL]*w_assets[j4, jjj, dy+1 , TL]
            counter += 1

            rowM[counter] = indnoh + (TL-5)*dy*(dy+1)*da + dy*dy*da + (jjj - 1)*da + j4
            colM[counter] = indret + polt[j4, jjj, dy+1, TL] + 1
            valM[counter] = sur[TL]*(1 - w_assets[j4, jjj, dy+1 , TL])
            counter += 1
        end  
    end
end


function transition_workers_3_1!(params, polt, w_assets, rowM, colM, valM, counter)
    """
    Transition matrix for working people who have not inherited so far and who do not inherit next period.
    Age 20 to 35 included
    
    Inputs:
    - params: dictionary of parameters
    - polt: policy function of workers (optimal asset index)
    - w_assets: weights assigned to two asset grid points on agrid surrounding the optimal savings

    Outputs- The following objects are updated:
    - rowM: add row indices of the nonzero elements of the transition matrix of workers
    - colM: add column indices of the nonzero elements of the transition matrix of workers
    - valM: add values of the nonzero elements of the transition matrix of workers
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of workers
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    polt[polt .>= da] .= da - 1

    for jj in 1:4 # loop over age
        for j0 in 1:dy # parents prod. state
            for jjj in 1:dy # agent's prod. state
                for j4 in 1:da # cash-on-hand
                    for j5 in 1:dy # agent's prod. state tmrw
                        rowM[counter] = (jj - 1)*dy*(dy)*da + (j0 - 1)*dy*da + (jjj - 1)*da + j4
                        colM[counter] = jj*dy*(dy)*da + (j0 - 1)*dy*da + (j5 - 1)*da + polt[j4, jjj, j0, jj]
                        valM[counter] = Qy[jjj, j5]*sur[jj + dage]*w_assets[j4, jjj, j0, jj]
                        counter += 1

                        rowM[counter] = (jj - 1)*dy*(dy)*da + (j0 - 1)*dy*da + (jjj - 1)*da + j4
                        colM[counter] = jj*dy*(dy)*da + (j0 - 1)*dy*da + (j5 - 1)*da + (polt[j4, jjj, j0, jj] + 1)
                        valM[counter] = Qy[jjj, j5]*sur[jj + dage]*(1 - w_assets[j4, jjj, j0, jj])
                        counter += 1
                    end
                end
            end
        end
    end
    
end


function transition_workers_3_2!(params, indnoh ,polt, w_assets, rowM, colM, valM, counter)
    """
    Transition matrix for working people who have not inherited so far and who do not inherit next period.
    Age 40 to 55 included
    
    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - polt: policy function of workers (optimal asset index)
    - w_assets: weights assigned to two asset grid points on agrid surrounding the optimal savings

    Outputs- The following objects are updated:
    - rowM: add row indices of the nonzero elements of the transition matrix of workers
    - colM: add column indices of the nonzero elements of the transition matrix of workers
    - valM: add values of the nonzero elements of the transition matrix of workers
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of workers
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    polt[polt .>= da] .= da - 1

    for jj in 5:TL-1 # loop over age
        for j0 in 1:dy # parents prod. state
            for jjj in 1:dy # agent's prod. state
                for j4 in 1:da # cash-on-hand
                    for j5 in 1:dy # agent's prod. state tmrw
                        rowM[counter] = indnoh + (jj - 5)*dy*(dy)*da + (j0 - 1)*dy*da + (jjj - 1)*da + j4
                        colM[counter] = indnoh + (jj - 4)*dy*(dy)*da + (j0 - 1)*dy*da + (j5 - 1)*da + polt[j4, jjj, j0, jj]
                        valM[counter] = Qy[jjj, j5]*sur[jj + dage]*w_assets[j4, jjj, j0, jj]
                        counter += 1

                        rowM[counter] = indnoh + (jj - 5)*dy*(dy)*da + (j0 - 1)*dy*da + (jjj - 1)*da + j4
                        colM[counter] = indnoh + (jj - 4)*dy*(dy)*da + (j0 - 1)*dy*da + (j5 - 1)*da + (polt[j4, jjj, j0, jj] + 1)
                        valM[counter] = Qy[jjj, j5]*sur[jj + dage]*(1 - w_assets[j4, jjj, j0, jj])
                        counter += 1
                    end
                end
            end
        end 
        
    end
end


function transition_workers_4_1!(params, indnoh, sopt, rowM, colM, valM, counter, beqdist)
    """
    Transition matrix for working people who have not inherited so far and will do so next period (only age 35)
    
    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - sopt: optimal asset choice of workers
    - beqdist: bequest distribution matrix

    Outputs- The following objects are updated:
    - rowM: add row indices of the nonzero elements of the transition matrix of workers
    - colM: add column indices of the nonzero elements of the transition matrix of workers
    - valM: add values of the nonzero elements of the transition matrix of workers
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of workers
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    for j0 in 1:dy # parents prod. state
        for jjj in 1:dy # agent's prod. state
            for j4 in 1:da # cash-on-hand
                polt_new = zeros(da)
                polt_new .= binary_search(agrid, (sopt[j4, jjj, j0, 4] .+ (agrid .- τ_b.*max.(agrid .- exb, zeros(size(agrid))))))
                polt_new[polt_new .>= da] .= da - 1

                wt_new = zeros(da)
                wt_new .= ((sopt[j4, jjj, j0, 4] .+ (agrid .- τ_b.*max.(agrid .- exb, zeros(size(agrid))))) .- agrid[polt_new])./(agrid[polt_new .+ 1] .- agrid[polt_new])
                for j5 in 1:dy # agent's prod. state tmrw
                    for j6 in 1:da # inherited assets
                        valM[counter + polt_new[j6] - 1] = valM[counter + polt_new[j6] - 1] +  sur[4]*(1 - sur[4 + dage])*Qy[jjj, j5]*beqdist[j6, j0, 1]*wt_new[j6]
                        valM[counter + polt_new[j6]] = valM[counter + polt_new[j6]] +  sur[4]*(1 - sur[4 + dage])*Qy[jjj, j5]*beqdist[j6, j0, 1]*(1 - wt_new[j6])
                    end
                    rowM[counter:counter+da-1] .= 3*dy*dy*da + (j0 - 1)*dy*da + (jjj - 1)*da + j4
                    colM[counter:counter+da-1] .= indnoh .+ dy*dy*da .+ (j5 - 1)*da .+ collect(1:da)
                    counter += da
                end
            end
        end
    end
    
end


function transition_workers_4_2!(params, indnoh, sopt, rowM, colM, valM, counter, beqdist)
    """
    Transition matrix for working people who have not inherited so far and will do so next period (age 40 to 55)
    
    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - sopt: optimal asset choice of workers
    - beqdist: bequest distribution matrix

    Outputs- The following objects are updated:
    - rowM: add row indices of the nonzero elements of the transition matrix of workers
    - colM: add column indices of the nonzero elements of the transition matrix of workers
    - valM: add values of the nonzero elements of the transition matrix of workers
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of workers
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens, 
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020, 
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    for jj in 5:TL-1 # loop over age
        tempint1 = indnoh + (jj-5)*dy*(dy+1)*da 
        tempint2 = indnoh + (jj-4)*dy*(dy+1)*da
        for j0 in 1:dy # loop over parents income state
            for jjj in 1:dy # loop over today's income
                for j4 in 1:da # loop over assets
                    polt_new = zeros(da)
                    polt_new .= binary_search(agrid, (sopt[j4, jjj, j0, jj] .+ (agrid .- τ_b.*max.(agrid .- exb, zeros(size(agrid))))))
                    polt_new[polt_new .>= da] .= da - 1

                    wt_new = zeros(da)
                    wt_new .= ((sopt[j4, jjj, j0, jj] .+ (agrid .- τ_b.*max.(agrid .- exb, zeros(size(agrid))))) .- agrid[polt_new])./(agrid[polt_new .+ 1] .- agrid[polt_new])

                    for j5 in 1:dy # loop over tomorrow's income
                        for j6 in 1:da # loop over inherited assets
                            valM[counter + polt_new[j6] - 1] = valM[counter + polt_new[j6] - 1] +  sur[jj]*(1 - sur[jj + dage])*Qy[jjj, j5]*beqdist[j6, j0, jj - 3]*wt_new[j6]
                            valM[counter + polt_new[j6]] = valM[counter + polt_new[j6]] +  sur[jj]*(1 - sur[jj + dage])*Qy[jjj, j5]*beqdist[j6, j0, jj - 3]*(1 - wt_new[j6])
                        end
                        rowM[counter:counter+da-1] .= tempint1 + (j0 - 1)*dy*da + (jjj - 1)*da + j4
                        colM[counter:counter+da-1] .= tempint2 + (j0 - 1)*dy*da + (j5 - 1)*da .+ collect(1:da)
                        counter += da 
                    end
                end
            end
        end 
    end
end

function transition_workers_4_3!(params, indnoh, indret, sopt, rowM, colM, valM, counter, beqdist)
    """
    Transition matrix for working people who have not inherited so far and will do so and retire next period (age 60)
    
    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - indret: total no. of states when agent was young and working
    - sopt: optimal asset choice of workers
    - beqdist: bequest distribution matrix

    Outputs- The following objects are updated:
    - rowM: add row indices of the nonzero elements of the transition matrix of workers
    - colM: add column indices of the nonzero elements of the transition matrix of workers
    - valM: add values of the nonzero elements of the transition matrix of workers
    - counter: counter to keep track of the number of nonzero elements of the transition matrix of workers
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens,
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020,
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    tempint1 = indnoh + (TL-5)*dy*(dy+1)*da

    for j0 in 1:dy # parents productivity state
        for jjj in 1:dy # income state
            for j4 in 1:da # asset state 
                polt_new = zeros(da)
                polt_new .= binary_search(agrid, (sopt[j4, jjj, j0, TL] .+ (agrid .- τ_b.*max.(agrid .- exb, zeros(size(agrid))))))
                polt_new[polt_new .>= da] .= da - 1

                wt_new = zeros(da)
                wt_new .= ((sopt[j4, jjj, j0, TL] .+ (agrid .- τ_b.*max.(agrid .- exb, zeros(size(agrid))))) .- agrid[polt_new])./(agrid[polt_new .+ 1] .- agrid[polt_new])

                for j6 in 1:da # inherited assets
                    valM[counter + polt_new[j6] - 1] = valM[counter + polt_new[j6] - 1] +  sur[TL]*beqdist[j6, j0, TL - 3]*wt_new[j6]
                    valM[counter + polt_new[j6]] = valM[counter + polt_new[j6]] +  sur[TL]*beqdist[j6, j0, TL - 3]*(1 - wt_new[j6]) 
                end

                rowM[counter:counter+da-1] .= tempint1 + (j0 - 1)*dy*da + (jjj - 1)*da + j4
                colM[counter:counter+da-1] .= indret .+ polt_new .+ 1
                counter += da
            end 
        end
        
    end

end

function beq_dist_update_1!(params, indnoh, rowM, colM, valM, stbeq)
    """
    Function to update the conditional characterstics of agents at age 40 conditional on their income.
    This function update the above using transition matrix.
    This will eventually become bequest distribution after we update it to age 60, and derive the marginal for asset distribution.

    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - rowM: row indices of the nonzero elements of the transition matrix of workers
    - colM: column indices of the nonzero elements of the transition matrix of workers
    - valM: values of the nonzero elements of the transition matrix of workers

    Outputs- The following objects are updated:
    - stbeq: updated bequest distribution matrix
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens,
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020,
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    for jj in dage:TL-1 # loop over age for which we need to update the joint dist.
        stbeq2 = similar(stbeq)
        stbeq2 .= 0.0

        # Case-1: First working people that have already inherited 
        for i in 4*da + (jj-5)*dy*dy*da*2 + 1:4*da*2 + (jj-4)*dy*dy*da*2
            tempint2 = rowM[i] - indnoh - (jj-5)*dy*(dy+1)*da
            tempint1 = colM[i] - indnoh - (jj-4)*dy*(dy+1)*da
            stbeq2[tempint1 ,: ] .= stbeq2[tempint1 ,: ] .+ valM[i].*stbeq[tempint2 ,: ]
        end

        # Case-2: People that have not inherited and keep not inheriting
        xx = 4*da + 4*dy*dy*da*2 + dy*da*2 + 4*dy*dy*da*dy*2
        for i in xx + (jj-5)*dy*dy*da*dy*2 + 1:xx + (jj-4)*dy*dy*da*dy*2
            tempint2 = rowM[i] - indnoh - (jj-5)*dy*(dy+1)*da
            tempint1 = colM[i] - indnoh - (jj-4)*dy*(dy+1)*da
            stbeq2[tempint1 ,: ] .= stbeq2[tempint1 ,: ] .+ valM[i].*stbeq[tempint2 ,: ] 
        end

        # Case-3: People that have not inherited and will inherit next period
        xx = 4*da + 4*dy*dy*da*2 + dy*da*2 + 8*dy*dy*da*dy*2 + dy*dy*da*dy*da
        for i in xx + (jj-5)*dy*dy*da*dy*da + 1:xx + (jj-4)*dy*dy*da*dy*da
            tempint2 = rowM[i] - indnoh - (jj-5)*dy*(dy+1)*da
            tempint1 = colM[i] - indnoh - (jj-4)*dy*(dy+1)*da
            stbeq2[tempint1 ,: ] .= stbeq2[tempint1 ,: ] .+ valM[i].*stbeq[tempint2 ,: ] 
        end

        stbeq .= stbeq2./sur[jj]
        
    end
    
end

function beq_dist_update_2!(params, indnoh, indret, nonzeroM, rowM, colM, valM, stbeq, stbeqret)
    """
    This function updates the conditional distribution of characterstics at age 60 conditional on income, for retirement period and later. 
    Inputs:
    - params: dictionary of parameters
    - indnoh: starting row index for people too young to have inherited
    - indret: total no. of states when agent was young and working
    - nonzeroM: number of nonzero elements in the transition matrix of workers
    - rowM: row indices of the nonzero elements of the transition matrix of workers
    - colM: column indices of the nonzero elements of the transition matrix of workers
    - valM: values of the nonzero elements of the transition matrix of workers
    - stbeq: updated bequest distribution matrix for workers (from last function)

    Outputs- The following objects are updated:
    - stbeqret: updated bequest distribution matrix for retirees
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens,
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020,
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params

    stbeqret2 = similar(stbeqret)
    stbeqret2 .= 0.0

    # Case-1: People that have already inherited and turn 65 next year (Still in the state-space of working people)
    for i in 4*da + 4*dy*dy*da*2 + 1: 4*da + 4*dy*dy*da*2 + dy*da # this last part does not have *2 coz in my case vfi of retirees does not have any interpolation coz optimal savings function lies on agrid by construction
        tempint2 = rowM[i] - indnoh - (TL-5)*dy*(dy+1)*da
        tempint1 = colM[i] - indret
        stbeqret2[tempint1 ,: ] .= stbeqret2[tempint1 ,: ] .+ valM[i].*stbeq[tempint2 ,: ] 
        
    end

    # Case-2: People that have not inherited and will inherit next period and retire next period
    xx = 4*da + 4*dy*dy*da*2 + dy*da*2 + 8*dy*dy*da*dy*2 + 5*dy*dy*da*dy*da
    for i in xx + 1:xx + nonzeroM
        tempint2 = rowM[i] - indnoh - (TL-5)*dy*(dy+1)*da
        tempint1 = colM[i] - indret
        stbeqret2[tempint1 ,: ] .= stbeqret2[tempint1 ,: ] .+ valM[i].*stbeq[tempint2 ,: ] 
    end
    stbeqret .= stbeqret2./sur[TL]

end

function beq_dist_update_3!(params, indret, nonzeroMT, rowM, colM, valM, rowMT, colMT, valMT, stbeqret, beqdist)
    """
    This function updates the conditional distribution of characterstics at age 60 conditional on income, for age 70 upto 90.
    Inputs:
    - params: dictionary of parameters
    - indret: total no. of states when agent was young and working
    - rowM: row indices of the nonzero elements of the transition matrix of workers
    - colM: column indices of the nonzero elements of the transition matrix of workers
    - valM: values of the nonzero elements of the transition matrix of workers
    - rowMT: row indices of the nonzero elements of the transition matrix of oldies (who are alive when 90)
    - colMT: column indices of the nonzero elements of the transition matrix of oldies (who are alive when 90)
    - valMT: values of the nonzero elements of the transition matrix of oldies (who are alive when 90)
    - stbeqret: updated bequest distribution matrix for retirees
    - nonzeroMT: number of nonzero elements in the transition matrix of oldies (who are alive when 90)

    Outputs- The following objects are updated:
    - stbeqret: updated bequest distribution matrix for retirees
    - beqdist: updated bequest distribution matrix for retirees (age 70 to 90)
    """
    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens,
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020,
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params


    # Update the bequest distribution for retirees: Age 70 to 85
    for jj in TL+1:T-1
        stbeqret2 = similar(stbeqret)
        stbeqret2 .= 0.0
        
        for i in (jj-TL-1)*da+1:(jj-TL)*da 
            tempint1 = colM[i] - indret - (jj-TL)*da
            tempint2 = rowM[i] - indret - (jj-TL-1)*da
            stbeqret2[tempint1 ,: ] .= stbeqret2[tempint1 ,: ] .+ valM[i].*stbeqret[tempint2 ,: ]
        end

        stbeqret .= stbeqret2./sur[jj]
        beqdist[ :, :, jj-dage-3] .= stbeqret
        
    end

    # Lastly we consider the saving dist. of agents who survived until their last period of life (Age 90)
    stbeqret2 = similar(stbeqret)
    stbeqret2 .= 0.0
    for i in 1:nonzeroMT 
        stbeqret2[colMT[i],:] = stbeqret2[colMT[i],:] .+ valMT[i].*stbeqret[rowMT[i],:]
    end
    beqdist[ :, :, T-dage-3] .= stbeqret2

end


function asset_by_age!(params, invm, indnoh, indret, nonzeroMT, rowMT, colMT, valMT, agew, agewT)
    """
    Asset distribution by age

    Inputs:
    - params: Dict, parameters of the model
    - invm: Array, invariant distribution of the state variables
    - indnoh: Int, starting row index for the people too young to have inherited
    - indret: Int, starting row index for the retirees
    - nonzeroMT: Int, number of nonzero elements in the transition matrix of oldies (who are alive when 90)
    - rowMT: Array, row indices of the nonzero elements of the transition matrix of oldies (who are alive when 90)
    - colMT: Array, column indices of the nonzero elements of the transition matrix of oldies (who are alive when 90)
    - valMT: Array, values of the nonzero elements of the transition matrix of oldies (who are alive when 90)

    Outputs: Objects that are modified in place
    - agew: Array, asset distribution by age for all ages before T
    - agewT: Array, asset distribution for age T
    """

    # Unpacking parameters
    @unpack tlength, dage, dy, def, TL, dsur, T, TR, β ,σ, ϕ₁, ϕ₂, γ, hbet, τₐ, τ_b, exb, ggovt, pens,
    maxa, mina, da, minc, maxc, dc, r, δ, gpop, σ_z, σ_h, y, Qy, Qyh, Qher, invher, inv40, jt4020,
    eff, sur, csur, numch, agrid, anet, cgrid, ϕ = params



    for j in 1:4 # age 20-35 
        for i in 1:dy*dy
            agew[ :,j] = agew[ :,j] .+ invm[(j-1)*dy*dy*da + (i-1)*da + 1:i*da + (j-1)*dy*dy*da]
        end
    end

    for j in 5:TL # age 40-60
        for i in 1:dy*(dy+1)
            agew[ :,j] = agew[ :,j] .+ invm[indnoh + (i-1)*da + (j-5)*dy*(dy+1)*da + 1: indnoh + i*da + (j-5)*dy*(dy+1)*da]
        end
    end

    for j in TL+1:T # AGE 65-90
        agew[ :,j] = invm[indret + (j-TL-1)*da + 1: indret + (j-TL-1)*da]
    end

    for i in 1:nonzeroMT 
        agewT[colMT[i]] = agewT[colMT[i]] + agew[rowMT(i),T]*valMT[i]
    end

end
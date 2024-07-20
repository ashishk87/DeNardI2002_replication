# Loading some required packages
using Distributions, Plots, PyPlot, LaTeXStrings, IterTools, QuantEcon, CSV, DataFrames, LinearAlgebra, Random, Statistics,Conda, CSVFiles
using Statistics, CSV, DataFrames, LinearAlgebra, Random,Conda, BenchmarkTools, MacroTools, Base.Threads
using Optim, GLM, FreqTables, ForwardDiff, PyCall, PlotlyJS, HTTP, CovarianceMatrices, StatsBase,Printf
using JuMP, Ipopt, NLopt, StatsBase, Econometrics, CategoricalArrays, PyCall,  RCall, PrettyTables,FixedEffectModels, XLSX, BlockArrays,BlockBandedMatrices, NLsolve, Dates, IJulia, Interact
using IJulia, Base.Threads, DelimitedFiles


function create_model(;
    tlength = 5, # 1) Period length
    dage = 5,    # 2) parent-child age differences
    dy = 4,      # 3) number of income states
    def = 45,    # 4) size of elflife
    TL = 9,      # 5) no. of working periods
    dsur = 65,   # 6)  alive from 20 to 65
    T = 14,      # 7)  maximim no. of periods alive

    # Preferences
    sig = 1.5, # 9) CRRA
    phi1 = -9.5, # 11) joy of bequueathing
    phi2 = 8.0, # 12) bequest parameter

    # Income and inheritance processes
    gam = 0.85,  # 13) income autoregre. coeff
    hbet = 0.677,  # 14) Inheritance autoregression
    pii = 3.1415926,  # 16) pi

    # Government & tech
    taua = 0.20, # 18) tax on capital income
    taub = 0.10, # 19) tax rate on estate
    exb = 8.5,  # 20) years avg earn for estate taxes exempt
    ggovt = 0.18, # 21) public expenditure
    pens = 0.40, # 22) pensions: repl rate*avg income

    # Grids
    maxa = 50.0, # 24) max level of capital
    mina = 0.0, # 25) min level of capital (do not change)
    da = 120, # 26) no. of grid points for capital
    # Policy functions will be dfined on a finer grid, then will interpolate
    dar = 120, # 27) no. of grid points for retired agents
    minc = 0.0, # 29) min cons proportion. do not change
    maxc = 1.0, # 30) max cons proportion. do not change
    
    # Polynomial extrapolation (for value fn)
    dagg = 30, # 31) no. of points added on the grid
    dpol =  2, # 32) dpol-1 is the degree of the polynomial

    # Relaxation parameter for the iteration on the govt b.c.
    relax = 0.7, # 33)

    # Kernel stuff
    numkernel = 31,  # 41) number of kernel functions
    numstdev = 0.5,   # 46)
    )

    # Function of parameters
    TR = T - TL # 8) retirement length
    bet = (0.9455)^tlength # 10) discount factor
    gpop = 1.012^tlength # 15) population growth rate
    r = (1.06)^(tlength - 1) # 17) exogenous interest rate
    delt = 1 - (94)^tlength # 23) depreciation rate, 6% per year
    # consumption grid for working agents (we are substituting a_{t+1})
    dc = dar # 28) no. of grid points for consumption

    # transition matrioces M and MT
    #=
     M:: number of non-zero elements in the transition matrix M
        20-25 yrs old: possible values today: (da*dy*dy)
        possible values tomorrow: (dy*2), (because of the interpolation and how we attribute the mass to 2 points)
        ! for 20 and 25 yrs old
        ! 30 yrs old: next period might inherit
        ! 35 to 55 yrs old
        ! 60 yrs old
        ! 65 to 85 yrs old
    =#
    nonzeroM = (da*dy*dy*3)*(dy*2)+(da*dy*dy)*(dy*2+dy*da)+  # 34)
               (da*dy*dy*4)*(dy*2+dy*da)+(da*dy*4)*(dy*2)+
               (da*dy*dy)*da+(da*dy)*2+(da*4)*2

    #=
    MT:: number of non-zero elements in the transition matrix MT
         ! the 90 yrs old guys
    =#
    nonzeroMT = da*2  # 35)

    # size of the invariant distribution invm
    ninvm=4*dy*dy*da+5*dy*(dy+1)*da+5*da  # 36)

    # starting row index -1 for the retired
    indret = 4*dy*dy*da+5*dy*(dy+1)*da  # 37)

    # starting row index -1 for the people too young to have inherited
    indnoh=dy*dy*da*4  # 38)

    # define dy*dy*da 
    dydyda=dy*dy*da  # 39)

    # define dy*(dy+1)*da 
    dydy1da=dy*(dy+1)*da  # 40)
    groupoint=da/(numkernel-1) # 41) number of points between

    #=
    ! the mean of 1 kernel and the mean of the next
    ! the first kernel is special (match probability at zero)
    =#
    middlepoint=(groupoint+1)/2  # 42)
    indmeank = zeros(Int64, numkernel-1) # 43)
    upint = zeros(Int64, numkernel-1) # 44)
    lowint = zeros(Int64, numkernel-1) # 45)
    stdevk = zeros(Float64, numkernel-1) # 47)
    meankernel = zeros(Float64, numkernel-1) # 48)
    kernelfun =  zeros(Float64, da, numkernel) # 49) value of the kernel
    expkernf  = zeros(Float64, numkernel*dy*(TL-3)) # 50)
    expkernf2 = zeros(Float64, numkernel*dy*(TL-3)) # 51)
    kercoeff  = zeros(Float64, numkernel, dy, TL-3) # 52)
    kercoeff2 = zeros(Float64, numkernel, dy, TL-3) # 53)
    kertemp =  zeros(Float64, numkernel, numkernel) # 54)
    kernelreg = zeros(Float64, numkernel, da) # 55)


    return (;
        tlength, dage, dy, def, TL, dsur, T, sig, phi1, phi2, gam, hbet, pii, taua, taub, exb, ggovt, pens, maxa, mina, da, dar, minc, maxc, dagg, dpol, relax, numkernel, numstdev, TR, bet, gpop, r, delt, dc, nonzeroM, nonzeroMT, ninvm, indret, indnoh, dydyda, dydy1da, groupoint, middlepoint, indmeank, upint, lowint, stdevk, meankernel, kernelfun, expkernf, expkernf2, kercoeff, kercoeff2, kertemp, kernelreg
    )
end


# The analogue of cumsumv in Julia is simply cumsum


function interplin(l, x, y, n, z, dpol, xinv)
    """
    ! interpolate, from values of the fn given as ordinates of
    ! input data points in an x-y plane and for a given set of
    ! x values (abscissas), the values of a single-valued fn y=f(x)
    ! input parameters are
    !     l  = # of data points (2 or greater)
    !     x  = array of dimension l storing the x values (abscissas)
    !          of input data points (in ascending order)
    !     y  = array of dimension l storing the y values (ordinates)
    !          of input data points
    !     n  = # of points at which interpolation of the y value (ordinate)
    !          (ordinate) is needed (1 or greater)
    !     z  = array of dimension n storing the x values (abscissas) of
    !             desired points (points that we need to interpolate over)
    !     dpol-1 = degree of polynomial extrapolation
    !     xinv = matrix for extrapolation computed in main code to save
    !            time
    ! output parameter
    !     v  = array of dimension n where the interpolated y values
    !          (ordinates) are to be displayed
    """
    # Ensure types are compatible with the Fortran example
    @assert l > 1
    @assert n > 0
    @assert dpol >= 1

    v = zeros(n)

    for i in 1:n
        diffs = abs.(z[i] .- x)
        ind = findmin(diffs)[2] # get index of the minimum difference. The [2] is to get the index since [1] is the minimum value
        diff = z[i] - x[ind]

        if abs(diff) < 1e-6
            v[i] = y[ind]
        else
            #= This block is executed if the interpolation point z(i) is beyond the highest x value in the dataset (x(l)) and is intended for extrapolation using a polynomial of degree dpol-1. The polynomial coefficients are presumably calculated elsewhere in the code (not shown in the snippet you provided), 
                and xinv is used along with these coefficients to extrapolate the value of y at z(i). This is a more complex operation than linear interpolation and is specific to cases where extrapolation is required beyond the provided dataset. =#
            if ind == l && diff >= 0
                zpoly = [z[i]^j for j in 0:(dpol-1)]
                v[i] = dot(zpoly, xinv * y[(end-dpol+1):end])
            else
                if diff < 0
                    v[i] = y[ind-1] + (z[i] - x[ind-1]) / (x[ind] - x[ind-1]) * (y[ind] - y[ind-1])
                else
                    v[i] = y[ind] + (z[i] - x[ind]) / (x[ind+1] - x[ind]) * (y[ind+1] - y[ind])
                end
            end
        end
    end

    return v
end


# The linspace function can be replaced by range frunction in Julia



"""
The main program begins here
"""

# Unpack the model parameters
tlength, dage, dy, def, TL, dsur, T, sig, phi1, phi2, gam, hbet, pii, 
taua, taub, exb, ggovt, pens, maxa, mina, da, dar, minc, maxc, dagg, dpol, relax, numkernel, 
numstdev, TR, bet, gpop, r, delt, dc, nonzeroM, nonzeroMT, ninvm, indret, indnoh, dydyda, dydy1da, groupoint, 
middlepoint, indmeank, upint, lowint, stdevk, meankernel, kernelfun, expkernf, expkernf2, kercoeff, kercoeff2, kertemp, 
kernelreg = create_model()

taul = 0.2823769  # 56) initialize tax rate on labor income
sigz  = sqrt(0.3) # 57) standard deviation of the income shock
sigh = sqrt(0.37)  # 58) standard deviation of the inheritance shock
scal = (sigz + sigh)/2  # 59) scale parameter for hermite quadrature

#=
! discretize income and inherit. process using HERMITE QUADRATURE
! want them on the same grid y, take some combi of the condit. distr.
! (normal with mean 0 and variance a linear combination of the two
! variances) define a scaling variable for the variance of the weigthing fn
!!INCOME PROCESS y, Qy
! load grid points and weights from a code that used the nag routines
! $markov.dat contains already sorted grid and wt data. write sort routine later$.
=#

# The following code is a direct translation of the Fortran code. But assumes 3 states for the income process. Uncomment this part if direct Fortran code translation is required

# In the paper 4 income states are used. So, I use that values directly.
# # Open the file
# filename = "markov.dat"
# # Read the contents
# data = readdlm(filename)

# # Assuming the file alternates between grid points and weights
# grid = data[1:3]  # Take every other element starting from the first
# wt = data[4:6]  # Take every other element starting from the second

# # Renormalize the weights so that it takes into account pi
# wt = wt/sqrt(2*pi)

# # Constructing the transition matrix for the markov process
# p = zeros(dy, dy) # Initializing
# for i in 1:dy
#     for j in 1:dy
#         p[i, j] = (scal / sigz) * exp(-((grid[j] - gam * grid[i])^2 * 
#                   scal^2) / (2 * sigz^2)) * wt[j] * exp(grid[j]^2 / 2)
#     end
#     p[i, :] .= p[i, :] ./ sum(p[i, :])  # Normalization for each row
# end

# productivity grid (or income grid)
y = [0.2594, 0.6513, 1.5355, 3.8547]
# Within lifetime transition matrix. Each element is the prob. of moving from income quartile i to income quartile j in 5 years coz that's the period length
Qy = [0.7132 0.2764 0.0104 0.0; 0.1467 0.6268 0.2210 0.0055; 0.0055 0.2210 0.6268 0.1467; 0.0 0.0104 0.2764 0.7132]
# inheritance transition matrix. Each element is the prob. of kids getting income quartile j from parents in income quartile i
Qyh = [0.5257 0.4228 0.0509 0.0007; 0.1461 0.5538 0.2825 0.0176; 0.0176 0.2825 0.5538 0.1461; 0.0007 0.0509 0.4228 0.5257]

#=
! matrix Qyh: inheritance of y: 40 years old parent-20 years old child
! want transmission  matrix Qher:
! 20 years old parent- her future 20 years old child
! (because want the invariant distribution at birth of these characterstics in population)
! have: c_20'=p_40'*Qyh, (i.e. children at 20 gets assigned what parents are at 40 times the inheritance transmission matrix) 
! We can rewrite the above as:
! c_20'=p_35'*Qy*Qyh, as
! and p_40'=p_35'*Qy. (i.e. parents at 40 gets assigned what parents are at 35 times the income transmission matrix)
! Keep doing it recursively until we get the following
=#
Qher = Qy * Qy * Qy * Qy * Qyh

# Compute the eigenvalues and eigenvectors of the transpose of Qher
eigenresult = eigen(transpose(Qher))

# Extract the principal eigenvalue and eigenvector
eigval = real(eigenresult.values[4])
eigvec = real(eigenresult.vectors[:, 4])

# Compute the invariant distribution for "her". distribution of characterstics at birth.
invher = eigvec / sum(eigvec)

# Distribution of the characteristics at 40 years of age: multiply Qy 4 times coz it gets from age 20 to 40 (20 ->25 ->30 ->35 ->40) and inverh is the invariant distribution of characterstics at age 20
inv40 = Qy * Qy * Qy * Qy * invher

# Calculate the joint distribution for parent's characterstics at age 40 and children's characterstics at age 20.
#=
So, (i, j)th entry of jt4020 is the distribution that the parent at age 40 is in the i-th income state and the child at age 20 is in the j-th income state.
=#
jt4020 = zeros(dy, dy)
for i in 1:dy
    for j in 1:dy
        jt4020[i, j] = Qyh[i, j] * inv40[i]
    end
end

# ! construct age-efficiency profile for workers
#Reading data from a file
eflife = readdlm("eflife.dat")
# Transform yearly age-efficiency profile according to our period length
eff = [sum(eflife[i:min(i+tlength-1, end)]) for i in 1:tlength:length(eflife)]
# Normalizing the transformed efficiency profile
eff *= TL / sum(eff)

# Construct conditional survival probabilities
# Read survival probabilities from file
surlife = readdlm("surlife.dat")
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


# ! statistics for the productivity markov process
# ! compute the gini index for productivity.
# ! the number of people for each income level is given by the invariant
# ! distribution of the markov process we are considering.

# invdy is the normalized largest eigenvector of Qy
eigen_result = eigen(transpose(Qy))
largest_index = argmax(abs.(eigen_result.values))
invdy = real(eigen_result.vectors[:, largest_index])/sum(real(eigen_result.vectors[:, largest_index]))
dist = zeros(dy+1)
dist[1] = 0.0
dist[2:end] = invdy
yy = zeros(dy+1)
yy[1] = 0.0
yy[2:end] = y
dyy = size(yy, 1)

# ! lorenz curve plot (cumulative productivity plotted against cumulative
# ! number of agents both normalized to one)
# cumsum of vectors
cumdist = cumsum(dist)
cumearn = cumsum(yy .* dist)
cumearn = cumearn / cumearn[end]

# ! compute the gini index
# ! the area under the lorenz curve is composed by trapeziums:
# !          |c
# ! b |      |
# !   |      |
# ! a |      |d
# ! whose area is (cd+ba)*ad/2.
# ! the formula for the gini index given by the ratio of the
# ! area between the 45 degree line and the the lorenz curve divided by
# ! the total area under the 45 degree line:
# ! (1/2 - \sum_{i=1}**{dyy-1}(e(i+1)+e(i))*(cumdist(i+1)-cumdist(i))/2)/(1/2)=
# ! 1-\sum_{i=1}**{dyy-1}(e(i+1)+e(i))*(cumdist(i+1)-cumdist(i))
area = 0.0
for i in 1:dyy-1
    area += (cumearn[i+1] + cumearn[i]) * (cumdist[i+1] - cumdist[i])
end
gini = 1 - area

# ! we also want to compute the gini index for earnings taking into account
# ! the age-productivity profile and how the inheritance process interacts
# ! with the earnings distribution. call it "adjusted gini". we thus need to
# ! perform some computation. form the transition matrix:
# !                                    age 1        age 2    .... age TL
# !                              prod1....prod dy  .....         ........
# ! age 1 productivity level 1
# !        productivity level 2
# !        ....................
# !        productivity level dy
# ! age 2 productivity level 1
# !        productivity level 2
# !        ....................
# !        productivity level dy
# ! ...........................
# ! age TL productivity level 1
# !        ....................
# !        productivity level dy
EM = zeros(Float64, TL*dy, TL*dy)
for i in 1:TL-1  # today's age
    for j in 1:dy  # today's productivity
        EM[(i-1)*dy+j, dy+(i-1)*dy+1:dy+i*dy] = sur[i] .* Qy[j, :]
    end
end
eq = zeros(TL*dy)

eigen_result = eigen(transpose(Qyh))
largest_index = argmax(abs.(eigen_result.values))
invher = real(eigen_result.vectors[:, largest_index])/sum(real(eigen_result.vectors[:, largest_index]))

eq[1:dy] = invher
# Now compute the invariant distribution
epsi = 1
eminv = copy(eq)
while epsi > 1e-8
    eminv1 = copy(eminv)
    eminv = ((transpose(EM) * eminv) / gpop)  + eq
    epsi = maximum(abs.(eminv1 - eminv))
end

# Renormalize the invariant distribution
eminv = eminv / sum(eminv)

# ! now the we have the invariant distribution over age and earnings
# ! we can compute the "adjusted" gini index that takes into account
# ! also the age efficiency profile.
# ! form a matrix with age on the rows and productivity on the columns
# ! containing the number of people in each cell in order to compute
# ! the age-efficiency profile IMPLIED by our income process (if we normalize
# ! the mean of the log income process wo do not have the right to other
# ! normalizations)
ageyn = transpose(reshape(eminv, (dy, TL)))
#! transform the matrix in terms of fractions of people of a given age
for i in 1:TL
    ageyn[i, :] = ageyn[i, :] / sum(ageyn[i, :])
end
#! compute the average age-productivity profile implied by our income process
endeff = zeros(TL)
endeff = ageyn*y
#! readjust the income process so that the implied age-efficiency profile
# ! is the one given in the eflife matrix
eff = eff./endeff
# ! adjusted gini index:construct y*eff for each possible
# ! y and eff and sort them in ascending order, keeping track of the
# ! corresponding number of people in that cell
aearn = zeros(TL*dy)
for i in 1:TL
    for j in 1:dy
        aearn[(i-1)*dy+j] = eff[i]*y[j]
    end
end
numsort = TL*dy
# !create vector of indexes to be permuted for the sorting fn
iperm = zeros(TL*dy)
iperm[1] = 1
for i in 2:numsort
    iperm[i] = iperm[i-1] + 1
end # simply a vector of 1:TL*dy

# Sort the aearn vector in ascending order
aearn1 = sort(aearn, rev=false)

aearn0 = zeros(TL*dy+1)
aearn0[1] = 0.0
aearn0[2:end] = aearn1

adist = zeros(TL*dy+1)
adist[1] = 0.0
adist[2:end] = eminv

# Compute the cumulative distribution
cumadist = cumsum(adist)  # cumulative no. of agents
cumaearn = cumsum(aearn0 .* adist / sum(aearn0 .* adist))  # cumulative adjusted productivity

# Compute the gini index
areaa = 0.0
for i in 1:TL*dy
    areaa += (cumaearn[i+1] + cumaearn[i]) * (cumadist[i+1] - cumadist[i])
end
agini = 1 - areaa  # This is wrong. Run the FORTRAN code to get the correct value

# ! compute now the gini index for productivity for just born people and
# ! for people who are TL years old
dist1 = zeros(dy+1)
dist1[1] = 0.0
dist1[2:end] = eminv[1:dy]

dist2 = zeros(dy+1)
dist2[1] = 0.0
dist2[2:end] = eminv[(TL-1)*dy+1:TL*dy]

yy = zeros(dy+1)
yy[1] = 0.0
yy[2:end] = y

# ! lorenz curve plot (cumulative productivity plotted against cumulative
# ! number of agents both normalized to one)
cumdist1 = cumsum(dist1/sum(dist1)) # cumulative distribution
cumdist2 = cumsum(dist2/sum(dist2)) # cumulative distribution
cumearn1 = cumsum(yy .* dist1) / sum(yy .* dist1) # cumulative productivity
cumearn2 = cumsum(yy .* dist2) / sum(yy .* dist2) # cumulative productivity

# ! compute the gini index for new agents
area = 0.0
for i in 1:dy
    area += (cumearn1[i+1] + cumearn1[i]) * (cumdist1[i+1] - cumdist1[i])
end
gini1 = 1 - area

# ! compute the gini index for Tl years old agents
area = 0.0
for i in 1:dy
    area += (cumearn2[i+1] + cumearn2[i]) * (cumdist2[i+1] - cumdist2[i])
end
gini2 = 1 - area

# ! compute per capita number of children
# ! we have: newborns_{t+1}=gpop*newborns_t (1); pop. law of motion
# ! newborns_t=numch*parents_t (2); numch is the per capita number of children
# ! (2): total number of children
# ! parents_t=newborns_{t-dage}*sur(1) (3)
# ! (3): the number of parents (at 25) equals the number of newborns 5 periods
# ! before times their probability of survival (one until they are 20 and sur(1)
# ! between 20 and 25).
# ! from (2) and (3): newborns_t=numch*newborns_{t-dage}*sur(1)  (4)
# ! from (1): newborns_t=gpop^(dage)*newborns_{t-dage}           (5)
# ! equate (4) and (5), simplify and get:
numch = gpop^dage/sur[1]
# !efc=1+numch*(/0.5,0.06,0.07,0.08,0.0,0.0,0.0,0.0/)
efc = zeros(TL-1)
efc .= 1.0

sur1 = copy(sur)
sur1[1] = 1.0
sur1[2:T] = sur[1:T-1]
pfrac = zeros(T)
for i in 1:T
    pfrac[i] = prod(sur1[1:i] ./ gpop)
end
pfrac = pfrac ./ sum(pfrac)

# ! aggregate labor income (which is exogenous).
aveinc = sum(pfrac[1:TL] .* eff .* endeff)
numret = sum(pfrac[TL+1:T])

# ! form assets grid for which the value fn will be defined
# ! how the grid is formed is crucial for later transformations
# ! to map one grid to another
# !a=real( (/ (i,i=0,da-1) /) )/real(da-1)*(maxa-mina)+mina
a = range(mina, sqrt(maxa), length=da)
a = a.^2
# ! define RECEIVED ASSETS NET OF ESTATE TAXES PER PERSON
# ! later on we will also divide it by csur(j0+4), j0 being age
anet = zeros(da)
anet[a .<= exb]  .= a[a .<= exb]./numch # numch divided coz assets are equally divided among children
anet[a .> exb] .= (a[a .> exb] .- taub .* (a[a .> exb] .- exb)) ./ numch


# ! kkkkkk KERNEL STUFF !kkkk**** modified apr23 2000
# ! find location of mean and each upper and lower bound
# ! ofreference intervals
numkernel = 31   # number of kernel functions
groupoint = da/(numkernel-1) # number of points between the mean of 1 kernel and the mean of the next
# The first kernel is special (match probability at zero)
middlepoint = trunc((groupoint+1)/2)
numstdev = 0.5

indmeank = zeros(Int64, numkernel-1)
upint = zeros(Int64, numkernel-1)
lowint = zeros(Int64, numkernel-1)
stdevk = zeros(Float64, numkernel-1)
for i in 1:numkernel-1
    indmeank[i] = middlepoint + (i-1)*groupoint
    upint[i] = groupoint + (i-1)*groupoint
    lowint[i] = 1 + (i-1)*groupoint
end
meankernel = a[indmeank]
stdevk = (a[upint] - a[lowint])/(2*numstdev)

# Initialize the kernel
kernelfun = zeros(da, numkernel)
kernelfun[1,1] = 1.0 # First kernel is degenerate on zero
# ! approximate the normal distributions that constitute the 
# ! kernel on our grid
dist = Normal(0, 1)
for j in 1:numkernel-1
    kernelfun[1, j+1] = ((a[2] - meankernel[j]) / (a[2] - a[1])) * 
                        (cdf(dist, (a[2] - meankernel[j]) / stdevk[j]) - 
                         cdf(dist, (a[1] - meankernel[j]) / stdevk[j])) +

                        (stdevk[j] / (sqrt(2.0*pi) * (a[2] - a[1]))) * 
                        (exp(-((a[2] - meankernel[j])^2 / (2.0 * stdevk[j]^2))) - 
                         exp(-((a[1] - meankernel[j])^2 / (2.0 * stdevk[j]^2))))
    
    for i in 2:da-1
        kernelfun[i, j+1] = ((meankernel[j] - a[i-1]) / (a[i] - a[i-1])) * 
                            (cdf(dist, (a[i] - meankernel[j]) / stdevk[j]) - 
                             cdf(dist, (a[i-1] - meankernel[j]) / stdevk[j])) -
                            (stdevk[j] / (sqrt(2.0*pi) * (a[i] - a[i-1]))) * 
                            (exp(-((a[i] - meankernel[j])^2 / (2.0 * stdevk[j]^2))) - 
                             exp(-((a[i-1] - meankernel[j])^2 / (2.0 * stdevk[j]^2)))) +
                            ((a[i+1] - meankernel[j]) / (a[i+1] - a[i])) * 
                            (cdf(dist, (a[i+1] - meankernel[j]) / stdevk[j]) - 
                             cdf(dist, (a[i] - meankernel[j]) / stdevk[j])) +
                            (stdevk[j] / (sqrt(2.0*pi) * (a[i+1] - a[i]))) * 
                            (exp(-((a[i+1] - meankernel[j])^2 / (2.0 * stdevk[j]^2))) - 
                             exp(-((a[i] - meankernel[j])^2 / (2.0 * stdevk[j]^2))))
    end
    
    kernelfun[da, j+1] = ((meankernel[j] - a[da-1]) / (a[da] - a[da-1])) * 
                         (cdf(dist,(a[da] - meankernel[j]) / stdevk[j]) - 
                          cdf(dist,(a[da-1] - meankernel[j]) / stdevk[j])) -
                         (stdevk[j] / (sqrt(2.0*pi) * (a[da] - a[da-1]))) * 
                         (exp(-((a[da] - meankernel[j])^2 / (2.0 * stdevk[j]^2))) - 
                          exp(-((a[da-1] - meankernel[j])^2 / (2.0 * stdevk[j]^2))))
end

# Think of this as variance-covariance matrix
kertemp = kernelfun' * kernelfun
# Create the inverse of the kernel matrix: inverse of variance co-variance matrix
kertemp = inv(kertemp)
# Create the kernel coefficients: weigh those kernels more which have lower variance
kenelreg = kertemp * kernelfun'

31*120   120*3  31*3
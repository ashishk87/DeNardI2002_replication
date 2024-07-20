# Loading some required packages
using Distributions, Plots, PyPlot, LaTeXStrings, IterTools, QuantEcon, CSV, DataFrames, LinearAlgebra, Random, Statistics,Conda, CSVFiles
using Statistics, CSV, DataFrames, LinearAlgebra, Random,Conda, BenchmarkTools, MacroTools, Base.Threads
using Optim, GLM, FreqTables, ForwardDiff, PyCall, PlotlyJS, HTTP, CovarianceMatrices, StatsBase,Printf
using JuMP, Ipopt, NLopt, StatsBase, Econometrics, CategoricalArrays, PyCall,  RCall, PrettyTables,FixedEffectModels, XLSX, BlockArrays,BlockBandedMatrices, NLsolve, Dates, IJulia, Interact
using IJulia, Base.Threads


mutable struct EconomicParameters
    # Parameters
    tlength::Int8 # Period length
    dage::Int8    # parent-child age differences
    dy::Int8      # number of income states
    def::Int8     # size of elflife
    TL::Int8      # no. of working periods
    dsur::Int8    # alive from 20 to 65
    T::Int8       # maximim no. of periods alive
    TR::Int8      # retirement length

    # Preferences
    sig::Float32 # CRRA
    bet::Float64 # discount factor
    phi1::Float32 # joy of bequueathing
    phi2::Float32 # bequest parameter

    # Income and inheritance processes
    gam::Float32  # income AR
    hbet::Float32 # Inheritance autoregression
    gpop::Float64 # population growth
    pi::Float64   # pi

    # Government & tech
    r::Float32  # exogenous interest rate
    taua::Float32 # tax on capital income
    taub::Float32 # tax rate on estate
    exb::Float32  # years avg earn for estate taxes exempt
    ggovt::Float32 # government spending
    pens::Float32 # pensions: repl rate*avg income
    delt::Float32 # depreciation rate, 6% per year

    # Grids
    maxa::Float32 # max level of capital
    mina::Float32 # min level of capital (do not change)
    da::Int8      # no. of grid points for capital
    dar::Int8     # no. of grid points for retired agents
    dc::Int8      # no. of grid points for consumption
    minc::Float32 # min cons proportion. do not change
    maxc::Float32 # max cons proportion. do not change

    # Polynomial extrapolation
    dagg::Int16  # no. of points added on the grid
    dpol::Int8   # dpol-1 is the degree of the polynomial

    # Relaxation parameter
    relax::Float32 # relaxation parameter for the iteration on the govt b.c.

    # New fields for transition matrices and related parameters
    nonzeroM::Int32
    nonzeroMT::Int16
    ninvm::Int32
    indret::Int32
    indnoh::Int16
    dydyda::Int16
    dydy1da::Int16

    # Kernel stuff
    numkernel::Int
    groupoint::Int
    middlepoint::Int
    indmeank::Vector{Int}
    upint::Vector{Int}
    lowint::Vector{Int}
    numstdev::Float32
    stdevk::Vector{Float32}
    meankernel::Vector{Float32}
    kernelfun::Matrix{Float32}
    expkernf::Vector{Float32}
    expkernf2::Vector{Float32}
    kercoeff::Array{Float32,3}
    kercoeff2::Array{Float32,3}
    kertemp::Matrix{Float32}
    kernelreg::Matrix{Float32}

    # Economic metrics and others
    mny::Float32
    vary::Float32
    mnyh::Float32
    varyh::Float32
    eigval::Float32
    epsi::Float32
    epsilon::Float32
    epsilon2::Float32
    area::Float32
    aarea::Float32
    gini::Float32
    agini::Float32
    gini1::Float32
    gini2::Float32
    numch::Float32
    aveinc::Float32
    numret::Float32
    taul::Float32
    taul2::Float32
    gbal::Float32
    netearn::Float32
    weal::Float32
    tempreal::Float32
    avgk::Float32
    beqtax::Float32
    wgini::Float32
    wginic::Float32
    fpoor::Float32
    trwealth::Float32
    gnp::Float32
    k2gnp::Float32
    kshare::Float32
    kotsum::Float32

    # age-cons. requirement profile
    efc::Vector{Float32}

    # Constructor
    function EconomicParameters()
        new(
            5, 5, 3, 45, 9, 65, 14, 5,  # Parameters
            1.5, 0.9455^5, -9.5, 8.0,   # Preferences
            0.83, 0.67, 1.012^5, 3.1415926, # Income & inheritance
            (1.06^5) - 1, 0.20, 0.10, 8.5, 0.18, 0.40, 1 - 0.94^5, # Gov & tech
            50.0, 0.0, 120, 120, 120, 0.0, 1.0,  # Grids
            30, 2,  # Polynomial extrapolation
            0.7 # Relaxation parameter


        )

        # Calculate dependent values
        dydyda = dy * dy * da
        dydy1da = dy * (dy + 1) * da
        nonzeroM = (da*dy*dy*3)*(dy*2) + (da*dy*dy)*(dy*2+dy*da) +
                (da*dy*dy*4)*(dy*2+dy*da) + (da*dy*4)*(dy*2) +
                (da*dy*dy)*da + (da*dy)*2 + (da*4)*2
        nonzeroMT = da * 2
        ninvm = 4*dy*dy*da + 5*dy*(dy+1)*da + 5*da
        indret = 4*dy*dy*da + 5*dy*(dy+1)*da
        indnoh = dy*dy*da * 4
        numkernel = 31
        groupoint = da ÷ (numkernel - 1)
        middlepoint = (groupoint + 1) ÷ 2
        # Initialize vectors/matrices with placeholder sizes
        indmeank = zeros(Int, numkernel-1)
        upint = zeros(Int, numkernel-1)
        lowint = zeros(Int, numkernel-1)
        stdevk = zeros(Float32, numkernel-1)
        meankernel = zeros(Float32, numkernel-1)
        kernelfun = zeros(Float32, da, numkernel)
        expkernf = zeros(Float32, numkernel*dy*(TL-3))
        expkernf2 = zeros(Float32, numkernel*dy*(TL-3))
        kercoeff = zeros(Float32, numkernel, dy, TL-3)
        kercoeff2 = zeros(Float32, numkernel, dy)
    end
end

params = EconomicParameters()

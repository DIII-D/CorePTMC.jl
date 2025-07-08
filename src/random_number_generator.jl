## -------- RANDOM GENERATOR ABSTRACT TYPE -------- ##

abstract type RandomGeneratorType{M<:Random.AbstractRNG} end # AbstractRNG is an abstract interface for Random Number Generators. This type provides a unified way to handle different kinds of RNGs, allowing them to be used interchangeably where an RNG is required, such as in functions for generating random numbers.
# Random.MersenneTwister (as Random.Xoshiro) are RNG subtypes based on MarsenneTwister and Xoshiro algorithms for Random Numbers Generation

RandomGenerator = AbstractScalarField{RandomGeneratorType}

## -------- RANDOM GENERATOR CONSTRUCTORS AND FUNCTIONS -------- ##

# constructor of instance of type RandomGeneratorType for scalars
AbstractScalarField{RandomGeneratorType}(Np::Int64) = ScalarField{Vector{Float64},RandomGeneratorType{MersenneTwister}}(zeros(Np))
# constructor of instance of type RandomGeneratorType for vectors
AbstractVector3D{RandomGeneratorType}(Np::Int64) = AbstractVector3D{RandomGeneratorType{MersenneTwister}}(zeros(Np), zeros(Np), zeros(Np))

# run RNG -> you get a number for each particle (vector of length Np)

function update!(r::ScalarField{Vector{Float64},T}, i::Int64) where {T<:RandomGeneratorType}  # updates i-th scalar component of an existing 1D vector of scalars with a random value 
    r.value[i] = rand()
end
function update!(r::ScalarField{Vector{Float64},T}) where {T<:RandomGeneratorType}            # updates an existing vector of scalars with a vector of random values
    rand!(r.value) # puts random value ∈ [0,1] in r.value
end

function random_vector(Np::Int64; seed=123)            # creates a 1D vector of random scalar values
    χ = AbstractScalarField{RandomGeneratorType}(Np)
    #m = MersenneTwister(seed)
    update!(χ)
    return χ
end

update!(r::ScalarField{Vector{Float64},<:RandomGeneratorType}, m; seed=123) = rand!(m, r.value) # scalar case

function update!(r::Vector3D{Vector{Float64},<:RandomGeneratorType}, m) # 3D vector case, substitutes vector components with m. 
    r.x = rand!(m, r.x)   # put value m in r.x, r.y and r.z
    r.y = rand!(m, r.y)
    r.z = rand!(m, r.z)
end

# E.g.: r = Vector3D{Vector{Float64}, RandomGeneratorType}(1); m = Xoshiro(seed); update!(r,m;seed) updates r with 3 random values starting from seed 123 and using Xoshiro RNG
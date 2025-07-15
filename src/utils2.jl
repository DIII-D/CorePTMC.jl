#=
Author: Jerome Guterl (guterlj@fusion.gat.com), Luca Cappelli (cappellil@fusion.gat.com)
Company: General Atomics
utils.jl (c) 2024=#


function suppress_output(f, args...)
redirect_stdout(devnull) do
    f(args...)
end
end


# Redefinition of Base.fill! function, now including also Vector3D as input and repeating Np times each component of a given vector  
Base.fill(B::Vector3D{T,U}, Np::Int64) where {T,U} = AbstractVector3D{U}(fill(B.x, Np::Int64), fill(B.y, Np::Int64), fill(B.z, Np::Int64))
function Base.fill!(B::Vector3D{Vector{Float64},U}, v::Vector{Float64}) where {U}
    fill!(B.x, v[1])
    fill!(B.y, v[2])
    fill!(B.z, v[3])
    nothing
end

# allow a scalar field f::ScalarField to be called as f[i] instead of f.value[i]
Base.getindex(s::ScalarField{Vector{Float64},U}, i::Int64) where {U} = s.value[i]
Base.fill(B::ScalarField{T,U}, Np::Int64) where {T,U} = AbstractScalarField{U}(fill(B.value, Np::Int64))

# tolist function
tolist(s::Symbol) = s == :all ? s : [s]
tolist(s::Vector) = s


cross_x(a, b, i) = a.y[i] * b.z[i] - a.z[i] * b.y[i]
cross_y(a, b, i) = -a.x[i] * b.z[i] + a.z[i] * b.x[i]
cross_z(a, b, i) = a.x[i] * b.y[i] - a.y[i] * b.x[i]

cross_x(a, b) = a.y * b.z - a.z * b.y
cross_y(a, b) = -a.x * b.z + a.z * b.x
cross_z(a, b) = a.x * b.y - a.y * b.x

# Custom function to calculate the Euclidean norm of our 3D vector structure
function vector_norm(v::Vector3D{Float64})
return sqrt(v.x^2 + v.y^2 + v.z^2)
end

function vector3D_mod(vec3D::Vector3D)
    return sqrt.(vec3D.x .^ 2 + vec3D.z .^ 2 + vec3D.z .^ 2)
end

dot(A::Vector3D, B::Vector3D, i::Int64) = B.x[i] * A.x[i] + A.y[i] * B.y[i] + A.z[i] * B.z[i]
norm(A::Vector3D, i::Int64) = sqrt(dot(A, A, i))

# Vector3D sum
function vector3D_sum(A::AbstractVector3D{U}, B::AbstractVector3D{U}) where U
x = A.x + B.x
y = A.y + B.y
z = A.z + B.z
return AbstractVector3D{U}(x,y,z)
end

function vector3D_sum(A::Vector{<:AbstractVector3D{U}}) where U
isempty(A) && throw(ArgumentError("Cannot sum an empty list of ScalarFields"))

T = typeof(A[1])  # Capture the concrete type
sumx_value = zero(A[1].x)
sumy_value = zero(A[1].y)
sumz_value = zero(A[1].z)

for vector3Dfield in A
    sumx_value += vector3Dfield.x
    sumy_value += vector3Dfield.y
    sumz_value += vector3Dfield.z
end

return T(sumx_value, sumy_value, sumz_value)
end


# ScalarField sum
function scalarfield_sum(A::AbstractScalarField{U}, B::AbstractScalarField{U}) where U
return AbstractScalarField{U}(A.value + B.value)
end

function scalarfield_sum(A::AbstractScalarField{U}, Bs::AbstractVector3D...) where U
# Check that all elements in Bs are AbstractScalarFields
if any(b -> !(b isa AbstractScalarField{U}), Bs)
    throw(ArgumentError("All Bs must be of type AbstractScalarField{U}"))
end

# Sum up the values from A and each element in Bs
sum_value = A.value
for b in Bs
    sum_value += b.value
end

return AbstractScalarField{U}(sum_value)
end

# vettore di vettori lo apro con ..., allora ho una lista di vettori
function scalarfield_sum(A::Vector{<:AbstractScalarField{U}}, Bs::Vector{<:AbstractScalarField{U}}...) where U
    isempty(A) && throw(ArgumentError("Cannot sum an empty list of ScalarFields"))

    # Capture the concrete type
    T = typeof(A[1])

    # Initialize the sum_value with the value of the first scalar field in A
    sum_result = [AbstractScalarField{U}(scalarfield.value) for scalarfield in A]

    # Sum the values of the scalar fields in A
    for (i, scalarfield) in enumerate(A)
    sum_result[i].value += scalarfield.value
    end

    # Iterate over the variadic arguments (Bs) and sum their values
    for b_vector in Bs
        for (i, scalarfield) in enumerate(b_vector)
        sum_result[i].value += scalarfield.value
        end
    end

    return sum_result
end


function scalarfield_sum(A::Vector{<:AbstractScalarField{U}}) where U
isempty(A) && throw(ArgumentError("Cannot sum an empty list of ScalarFields"))

T = typeof(A[1])  # Capture the concrete type
sum_value = zero(A[1].value)

for scalarfield in A
    sum_value += scalarfield.value
end

return T(sum_value)
end


function hcat_vector3D(A::AbstractVector3D{U}, Bs::AbstractVector3D{U}...) where U
# Ensure all inputs have matching types
for B in Bs
    @assert typeof(A.x) == typeof(B.x) "typeof(A.x) must be equal to typeof(B.x)"
    @assert typeof(A.y) == typeof(B.y) "typeof(A.y) must be equal to typeof(B.y)"
    @assert typeof(A.z) == typeof(B.z) "typeof(A.z) must be equal to typeof(B.z)"
end

# Concatenate x components
if A.x isa Vector{<:Vector}
    N = length(A.x)
    x = [hcat(A.x[i], (B.x[i] for B in Bs)...) for i in 1:N]
else
    x = hcat(A.x, (B.x for B in Bs)...)
end

# Concatenate y components
if A.y isa Vector{<:Vector}
    N = length(A.y)
    y = [hcat(A.y[i], (B.y[i] for B in Bs)...) for i in 1:N]
else
    y = hcat(A.y, (B.y for B in Bs)...)
end

# Concatenate z components
if A.z isa Vector{<:Vector}
    N = length(A.z)
    z = [hcat(A.z[i], (B.z[i] for B in Bs)...) for i in 1:N]
else
    z = hcat(A.z, (B.z for B in Bs)...)
end

return AbstractVector3D{U}(x, y, z)
end

function vcat_vector3D(A::AbstractVector3D{U}, Bs::AbstractVector3D{U}...) where U
# Ensure all inputs have matching types
for B in Bs
    @assert typeof(A.x) == typeof(B.x) "typeof(A.x) must be equal to typeof(B.x)"
    @assert typeof(A.y) == typeof(B.y) "typeof(A.y) must be equal to typeof(B.y)"
    @assert typeof(A.z) == typeof(B.z) "typeof(A.z) must be equal to typeof(B.z)"
end

# Concatenate x components
if A.x isa Vector{<:Vector}
    N = length(A.x)
    x = [vcat(A.x[i], (B.x[i] for B in Bs)...) for i in 1:N]
else
    x = vcat(A.x, (B.x for B in Bs)...)
end

# Concatenate y components
if A.y isa Vector{<:Vector}
    N = length(A.y)
    y = [vcat(A.y[i], (B.y[i] for B in Bs)...) for i in 1:N]
else
    y = vcat(A.y, (B.y for B in Bs)...)
end

# Concatenate z components
if A.z isa Vector{<:Vector}
    N = length(A.z)
    z = [vcat(A.z[i], (B.z[i] for B in Bs)...) for i in 1:N]
else
    z = vcat(A.z, (B.z for B in Bs)...)
end

return AbstractVector3D{U}(x, y, z)
end


function vcat_vector3D(vec::Vector{<:AbstractVector3D{U}}) where U
@assert !isempty(vec) "Input vector of AbstractVector3D objects is empty."
first_obj = vec[1]

# Process the x-field.
if first_obj.x isa Vector{<:Vector}
    # Assume all objects have the same number of sub-vectors.
    N = length(first_obj.x)
    # For each index, concatenate the corresponding sub-vector from each object.
    x = [reduce(vcat, [obj.x[i] for obj in vec]) for i in 1:N]
else
    # Otherwise, simply concatenate the x fields from all objects.
    x = vcat((obj.x for obj in vec)...)
end

# Process the y-field.
if first_obj.y isa Vector{<:Vector}
    N = length(first_obj.y)
    y = [reduce(vcat, [obj.y[i] for obj in vec]) for i in 1:N]
else
    y = vcat((obj.y for obj in vec)...)
end

# Process the z-field.
if first_obj.z isa Vector{<:Vector}
    N = length(first_obj.z)
    z = [reduce(vcat, [obj.z[i] for obj in vec]) for i in 1:N]
else
    z = vcat((obj.z for obj in vec)...)
end

return AbstractVector3D{U}(x, y, z)
end



function vcat_scalarfield(A::AbstractScalarField{U}, Bs::AbstractScalarField{U}...) where U
    value = vcat(A.value, (B.value for B in Bs)...)
    return AbstractScalarField{U}(value)
end

function vcat_scalarfield(vec::Vector{<:AbstractScalarField{U}}) where U
    @assert !isempty(vec) "Input vector of AbstractScalarField objects is empty."
    value = vcat((obj.value for obj in vec)...)
    return AbstractScalarField{U}(value)
end



vcat_preserve_type(A::T, B::T) where {T} = T(vcat(A,B))

function vcat_preserve_type(A::Weights{T}, B::Weights{T}) where {T}  # special dispatch for weights
    return Weights(vcat(A, B))
end

# Generic version for any type T (assuming T is callable and can be constructed from a vector)
function vcat_preserve_type(vec::Vector{T}) where T
    return T(vcat(vec...))
end

# Specialized version for Weights (if needed)
function vcat_preserve_type(vec::Vector{<:Weights{T}}) where T
    return Weights(vcat(vec...))
end

function vcat_selector(A, B)
    if A isa AbstractVector3D && B isa AbstractVector3D
        merged = vcat_vector3D(A, B)
    elseif A isa AbstractScalarField && B isa AbstractScalarField
        merged = vcat_scalarfield(A, B)
    elseif A isa RandomGenerators && B isa RandomGenerators
        merged = A
    else
        merged = vcat_preserve_type(A, B)
    end
end

function vcat_selector(A, Bs...)
    if A isa AbstractVector3D && all(B -> B isa AbstractVector3D, Bs)
        merged = vcat_vector3D(A, Bs...)
    elseif A isa AbstractScalarField && !(A isa ScalarField{<:Array, SpatialParticleCountType})
        merged = vcat_scalarfield(A, Bs...)
    elseif A isa RandomGenerators && all(B -> B isa RandomGenerators, Bs)
        merged = A  # Assuming merging multiple RandomGenerators results in returning the first one
    elseif A isa ScalarField{<:Array, SpatialParticleCountType}
        merged = scalarfield_sum_sum(A, Bs...)
    else
        merged = vcat_preserve_type(A, Bs...)
    end
    return merged
end

function vcat_selector(A::Vector)
    println("A::Vector")
    if A isa Vector{<:AbstractVector3D} 
        merged = vcat_vector3D(A)
    elseif A isa Vector{<:AbstractScalarField} && !(A isa Vector{<:ScalarField{<:Array, SpatialParticleCountType}})
        merged = vcat_scalarfield(A)
    elseif A isa Vector{<:RandomGenerators}
        merged = A[1]
    elseif A isa Vector{<:ScalarField{<:Array, SpatialParticleCountType}}
        println("A isa SpatialParticleCountType")
        merged = scalarfield_sum(A)
    else
        merged = vcat_preserve_type(A)
    end
    return merged
end

function vcat_selector(A::Vector, Bs...)
    if A isa Vector{<:AbstractVector3D} && all(B -> B isa Vector{<:AbstractVector3D}, Bs)
        merged = vcat_vector3D.(A, Bs...)
    elseif A isa Vector{<:AbstractScalarField} && all(B -> B isa Vector{<:AbstractScalarField}, Bs)
        merged = vcat_scalarfield.(A, Bs...)
    elseif A isa Vector{<:RandomGenerators} && all(B -> B isa Vector{<:RandomGenerators}, Bs)
        merged = A[1]  # Assuming we always take the first one
    elseif A isa Vector{<:Vector} && all(B -> B isa Vector{<:Vector}, Bs)
        merged = vcat.(A, Bs...)
    else
        merged = vcat_preserve_type(A, Bs...)
    end
    return merged
end



function merge_structures(structA, structB) # merge all fields of 2 structures of the same type

    @assert typeof(structA) == typeof(structB) "structA and structB must be of the same type"

    field_names = fieldnames(typeof(structA))

    merged_fields = map(fn -> vcat_selector(getfield(structA, fn), getfield(structB, fn)), field_names)

    return typeof(structA)(merged_fields...)  # returns merged structure
end

## Integration functions

function trapz(x, y)
    n = length(y)
    result = 0

    for i in 2:n
        result += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
    end

    return result
end

#=
function trapz(x, y; dims=1)
    # Check if y is 1D
    if ndims(y) == 1
        n = length(y)
        result = 0
        for i in 2:n
            result += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
        end
        return result
    end

    # For multidimensional arrays
    if dims > ndims(y)
        throw(ArgumentError("dims cannot exceed the number of dimensions of y"))
    end

    # Permute the specified dimension to the first dimension
    permuted_dims = [dims; setdiff(1:ndims(y), dims)]
    y_permuted = permutedims(y, permuted_dims)

    # Integrate along the first dimension
    size_y = size(y_permuted)
    n = size_y[1]
    result = zeros(size_y[2:end]...)  # Result array without the first dimension
    for i in 2:n
        slice = 0.5 * (y_permuted[i, :, :] .+ y_permuted[i - 1, :, :]) .* (x[i] - x[i - 1])
        result .= result .+ slice
    end

    # Reverse the permutation
    inverse_permutation = invperm(permuted_dims)
    return permutedims(result, inverse_permutation)
end
=#

function cumtrapz(x, y)
    n = length(y)
    result = zeros(Float64, n)

    for i in 2:n
        result[i] = result[i - 1] + 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
    end

    return result
end

# export functions from LPTMC to Julia's Main
macro exportAll()
    inI = []
    out = []
    res = []

    for name in Base.names(__module__, all=true)
        push!(out, :(export $(name)))
    end

    for name in Base.names(__module__, all=false, imported=true)
        push!(inI, :(export $(name)))
    end

    for name in Base.names(Main.Base, all=false, imported=true)
        push!(inI, :(export $(name)))
    end

    for name in Base.names(Main.Core, all=false, imported=true)
        push!(inI, :(export $(name)))
    end

    for e in out
        if e in inI || e == :(export include) || '#' in string(e)
        continue
        end
        push!(res, :($(e)))
    end
    #= Do not export things that we already have =#
    ret = Expr(:block, res...)
    return ret
end
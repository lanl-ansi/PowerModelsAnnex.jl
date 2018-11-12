# This will soon become its own package, but is not yet.
module Units

using Unitful
using Missings

import Compat: @__MODULE__
using Compat.Dates

import Compat.Statistics: mean

export round, mean, asqtype, ustrip, fustrip, UnitfulMissing

import Base: round, *, convert
import Unitful: @derived_dimension, ustrip, uconvert, Quantity, W, hr, J, 𝐋, 𝐌, 𝐓

UnitfulMissing = Union{<:Unitful.Quantity, Missings.Missing}

@derived_dimension PowerHour 𝐋^2*𝐌*𝐓^-2
@unit Wh "Wh" WattHour 3600J true

@derived_dimension ReactivePowerHour 𝐋^2*𝐌*𝐓^-2
@unit VARh "VARh" VARHour 3600J true

@dimension Money "Money" Currency
@refunit USD "USD" Currency Money false

round(x::AbstractArray{<:Quantity}, r::Int) = round(ustrip(x), r) * unit(eltype(x))
round(x::T, r::Int) where T <: Quantity = round(ustrip(x), r) * unit(x)
mean(x::AbstractArray{<:Quantity}, r::Int) = mean(ustrip(x), r) * unit(eltype(x))
asqtype(x::T) where T <: Unitful.Units = typeof(1.0*x)
#ustrip{T<:Unitful.Quantity}(x::DataArrays.DataArray{T}) = map(t -> ustrip(t), x)
ustrip(x::Missings.Missing) = x
fustrip(x::Array{T}) where T = map(t -> ustrip(t), x)

Unitful.register(@__MODULE__)

*(y::Missings.Missing, x::T) where T <: Unitful.FreeUnits = y

*(x::Unitful.Units, y::Number) = *(y, x)

# Handle working with `Period`s
*(x::Unitful.Units, y::Period) = *(y, x)
*(x::Period, y::Unitful.Units) = convert(y, x)
function convert(a::Unitful.Units, x::Period)
    sec = Dates.value(Dates.Second(x))
    uconvert(a, (sec)u"s")
end

# Methods to drop
# Exist to test that (offsets)u"hr" should work the same way
dt2umin(t::AbstractArray{Dates.Minute}) = Dates.value.(t).*u"minute"

# Some gymnastics needed to get this to work at run-time.
# Sourced from https://github.com/ajkeller34/UnitfulUS.jl
const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(Units)
end

end

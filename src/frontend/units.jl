# This will soon become its own package, but is not yet.
__precompile__(false)
module Units

using Compat.Dates
using Missings

import Base: round, *, convert
import Compat.Statistics: mean

import Unitful
import Unitful: ustrip, J, W, hr, 𝐋, 𝐌, 𝐓
using Unitful: @unit, @dimension, @derived_dimension, @refunit, @u_str, uconvert, Quantity

export round, mean, asqtype, ustrip, fustrip, UnitfulMissing

const UnitfulMissing = Union{<:Unitful.Quantity, Missings.Missing}

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

# Re enable pre-compilation when this becomes its own package
# Some gymnastics needed to get this to work at run-time.
# Sourced from https://github.com/ajkeller34/UnitfulUS.jl
# const localunits = Unitful.basefactors
# const localpromotion = Unitful.promotion
# function __init__()
#     merge!(Unitful.basefactors, localunits)
#     merge!(Unitful.promotion, localpromotion) # only if you've used @dimension
#     Unitful.register(Units)
# end

end

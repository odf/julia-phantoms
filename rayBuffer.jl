abstract type Buffer{T <: Real} end


function update!(buffer::Buffer{T}, near::T, far::T, value::T) where T
    if length(buffer.empty) == 0
        push!(buffer.empty, (-Inf, Inf))
    end

    if near >= far
        return
    end

    i = 1
    while i <= length(buffer.empty)
        (enear, efar) = buffer.empty[i]
        if near <= enear < efar <= far
            updateValue!(buffer, enear, efar, value)
            splice!(buffer.empty, i)
            i -= 1
        elseif enear < near < far < efar
            updateValue!(buffer, near, far, value)
            buffer.empty[i] = (enear, near)
            insert!(buffer.empty, i + 1, (far, efar))
            i += 1
        elseif enear < near < efar <= far
            updateValue!(buffer, near, efar, value)
            buffer.empty[i] = (enear, near)
        elseif near <= enear < far < efar
            updateValue!(buffer, enear, far, value)
            buffer.empty[i] = (far, efar)
        end
        i += 1
    end
end


mutable struct ValueBuffer{T} <: Buffer{T}
    value::T
    empty::Vector{Tuple{T, T}}
    ValueBuffer{T}() where T = new(0, Vector{Tuple{T, T}}())
end


function clear!(buffer::ValueBuffer{T}) where T
    buffer.value = 0
    resize!(buffer.empty, 0)
end


function updateValue!(buffer::ValueBuffer{T}, near::T, far::T, value::T) where T
    buffer.value += value * (far - near)
end


function get(buffer::ValueBuffer{T}) where T
    return buffer.value
end


mutable struct AdditiveValueBuffer{T} <: Buffer{T}
    value::T
    AdditiveValueBuffer{T}() where T = new(0)
end


function clear!(buffer::AdditiveValueBuffer{T}) where T
    buffer.value = 0
end


function update!(
    buffer::AdditiveValueBuffer{T}, near::T, far::T, value::T
) where T
    buffer.value += value * (far - near)
end


function get(buffer::AdditiveValueBuffer{T}) where T
    return buffer.value
end


mutable struct IntervalBuffer{T} <: Buffer{T}
    intervals::Vector{Tuple{T, T, T}}
    empty::Vector{Tuple{T, T}}
    IntervalBuffer{T}() where T = new(Vector{Tuple{T, T, T}}(),
                                      Vector{Tuple{T, T}}())
end


function clear!(buffer::IntervalBuffer{T}) where T
    resize!(buffer.intervals, 0)
    resize!(buffer.empty, 0)
end


function updateValue!(
    buffer::IntervalBuffer{T}, near::T, far::T, value::T
) where T
    push!(buffer.intervals, (near, far, value))
end


function get(buffer::IntervalBuffer{T}) where T
    return copy(buffer.intervals)
end


mutable struct AdditiveIntervalBuffer{T} <: Buffer{T}
    intervals::Vector{Tuple{T, T, T}}
    AdditiveIntervalBuffer{T}() where T = new(Vector{Tuple{T, T, T}}())
end


function clear!(buffer::AdditiveIntervalBuffer{T}) where T
    resize!(buffer.intervals, 0)
end


function update!(
    buffer::AdditiveIntervalBuffer{T}, near::T, far::T, value::T
) where T
    push!(buffer.intervals, (near, far, value))
end


function get(buffer::AdditiveIntervalBuffer{T}) where T
    return copy(buffer.intervals)
end

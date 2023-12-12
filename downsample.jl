function linearDownsamplingKernel(n::Int)
    t = [n - abs(i - n) for i in 1 : 2 * n - 1] / n^2

    return [x * y for x in t, y in t]
end


function downsample(img::Array{T, 2}, factor::Int)::Array{T, 2} where T <: Real
    n, m = [div(k - 1, factor) for k in size(img)]
    out = Array{T}(undef, n, m)

    kernel = linearDownsamplingKernel(factor)

    for x in 1 : n
        xbase = (x - 1) * factor
        for y in 1 : m
            ybase = (y - 1) * factor

            t = zero(T)
            for i in 1 : 2 * factor - 1
                for j in 1 : 2 * factor - 1
                    t += kernel[i, j] * img[xbase + i, ybase + j]
                end
            end
            out[x, y] = t
        end
    end

    return out
end

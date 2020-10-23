"""
Basic methods for solving `f(x) = 0` equations and finding `argmin f(x)`
    
See [`bisection`](@ref) and [`goldensection`](@ref).
"""
module Dichotomy

export dichotomy, bisection, goldensection

"""
    bisection(f, (a, b); [tol=eps()])
Finds root of equation `f(x) = 0` on `a < x < b`.

Method works only if `f(a) * f(b) < 0`.

## Examples

```julia-repl
julia> bisection(sin, (2, 4))
3.141592653589793
```
"""
bisection(f, interval; args...) =
    bisection(f, extrema(interval)...; args...)

function bisection(f, a::Real, b::Real; tol::Real = eps())
    sa = sign(f(a))
    sb = sign(f(b))

    if sa * sb ≥ 0
        throw(ArgumentError("f(a) * f(b) ≥ 0"))
    end

    if sa > 0
        a, b = b, a
    end

    maxiter = ceil(Int, log(2, abs(b-a) / abs(tol)))

    for _ in 1:maxiter
        c = a / 2 + b / 2

        if f(c) < 0
            a = c
        else
            b = c
        end
    end

    return a / 2 + b / 2
end

const dichotomy = bisection

"""
    goldensection(f, (a, b); [tol=eps()])
Finds argmin of function `f(x)` on `a < x < b`.

## Examples

```julia-repl
julia> goldensection(x-> (x-2)^2, (-4, 4))
2.0
```
"""
goldensection(f, interval; args...) =
    goldensection(f, extrema(interval)...; args...)

function goldensection(f, a::Real, b::Real; tol::Real = eps())
    r = (3 - √5) / 2
    s = 1 - r

    α = s * a + r * b
    fα = f(α)

    maxiter = ceil(Int, log(s, abs(tol) / abs(b-a)))

    for _ in 1:maxiter
        β = s * b + r * a
        fβ = f(β)

        if fα < fβ
            b = a
            a = β
        else
            a = α
            α = β
            fα = fβ
        end
    end

    return α
end

end # module

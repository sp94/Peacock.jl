using Setfield
using Base.Iterators


"""
    AbstractShape

Supertype for shapes.

Subtypes shouuld have the following properties:
  ep::ComplexF64
  mu::ComplexF64
  transform::ShapeTransform

Subtypes should override the following methods:
  is_inside(shape, x, y)

Shapes can be easily translated and rotated, for example
```
shape |> translate(1,0) |> rotate(60)
```
will translate a shape by `(dx,dy)=(1,0)` and then rotate it by 60 degrees.
"""
abstract type AbstractShape end

# Implement Base.iterate(shape) so we can flatten a nested list of shapes
import Base.iterate
function Base.iterate(shape::AbstractShape)
    return shape, nothing
end
function Base.iterate(shape::AbstractShape, state)
    return nothing
end


"""
    sample(shapes::Array{<:AbstractShape,1}, x, y, fieldname::Symbol)

Return the value of the requested field (`ep` or `mu`) at `(x,y)`, working
in reverse order such that the last shape is on top.

If there is no `Background`, assume, `ep=mu=1`, for example,
```
shapes = [Circle(ep=4)]  # assumes background is ep=1, mu=1
shapes = [Background(ep=2), Circle(ep=4)]  # overrides default background
```
"""
function sample(shapes::AbstractArray{<:AbstractShape,1}, x, y, fieldname::Symbol)
    for shape in reverse(shapes)
        if is_inside(shape, x, y)
            return getfield(shape, fieldname)
        end
    end
    return 1.0 # default background
end


"""
    is_inside(shape::AbstractShape, x::Real, y::Real)

Returns whether `(x,y)` is inside the shape, taking
into account the `transform` of the shape.
"""
function is_inside(shape::AbstractShape)
    throw("Subtypes of AbstractShape should override is_inside(shape,x,y).")
end


"""
    ShapeTransform(origin::Vector{Float64}, axes::Matrix{Float64})

Stores the transform of a shape.

The `axes` matrix stores the x axis and y axis of the transformed
shape as columns of the matrix, ie `axes = [x_axis yaxis]`.
"""
struct ShapeTransform
    origin::Vector{Float64}
    axes::Matrix{Float64} # = [x_axis y_axis]
end


function ShapeTransform(origin::Vector{<:Real}, rotation_degrees::Real)
    axes = rotation_matrix(rotation_degrees)
    return ShapeTransform(origin, axes)
end


import Base.\
function \(t::ShapeTransform, xy::Vector{<:Real})
    return t.axes \ (xy - t.origin)
end


"""
    translate(transform::ShapeTransform, dx::Real, dy::Real)

Translate the origin of the `transform` by `(dx,dy)`.
"""
function translate(transform::ShapeTransform, dx::Real, dy::Real)
    return ShapeTransform(transform.origin+[dx;dy], transform.axes)
end


"""
TODO
"""
function translate(shape::AbstractShape, dx, dy)
    new_transform = translate(shape.transform, dx, dy)
    return @set shape.transform = new_transform
end


"""
    translate(dx, dy)

An alias for x->translate(x, dx, dy).

This is useful for piping shapes into transforms, for example,
```
Circle(radius) |> translate(1,0)
```
"""
function translate(dx, dy)
    return x->translate(x, dx, dy)
end


"""
TODO
"""
function rotate(transform::ShapeTransform, num_degrees::Real, origin::Union{Vector{<:Real},Nothing}=nothing)
    if origin == nothing
        origin = transform.origin
    end
    R = rotation_matrix(num_degrees)
    new_origin = R*(transform.origin-origin) + origin
    new_axes = R*transform.axes
    return ShapeTransform(new_origin, new_axes)
end


"""
TODO
"""
function rotate(shape::AbstractShape, num_degrees::Real, origin::Union{Vector{<:Real},Nothing}=nothing)
    new_transform = rotate(shape.transform, num_degrees, origin)
    return @set shape.transform = new_transform
end


"""
    rotate(num_degrees, origin=shape.origin)

An alias for x->translate(x, dx, dy).

For example,
```
Ellipse(radius1,radius2) |> rotate(90)
```
"""
function rotate(num_degrees::Real, origin::Union{Vector{<:Real},Nothing}=nothing)
    return x->rotate(x, num_degrees, origin)
end


"""
    Background(ep::ComplexF64=1, mu::ComplexF64=1)

Create a `Background` with permittivity `ep` and permeability `mu`.
"""
struct Background <: AbstractShape
    ep::ComplexF64
    mu::ComplexF64
    function Background(ep=1, mu=1)
        return new(ep, mu)
    end
end

function is_inside(background::Background, x::Real, y::Real)
    return true
end


"""
TODO
"""
struct Circle <: AbstractShape
    ep::ComplexF64
    mu::ComplexF64
    radius::Float64
    transform::ShapeTransform
end


"""
TODO
"""
function Circle(radius; ep=1, mu=1, origin=[0,0], angle=0)
    return Circle(ep, mu, radius, ShapeTransform(origin,angle))
end

function is_inside(circle::Circle, x::Real, y::Real)
    x, y = circle.transform \ [x; y]
    return x^2 + y^2 <= circle.radius^2
end


"""
TODO
"""
struct Ellipse <: AbstractShape
    ep::ComplexF64
    mu::ComplexF64
    radius1::Float64
    radius2::Float64
    transform::ShapeTransform
end


"""
TODO
"""
function Ellipse(radius1, radius2; ep=1, mu=1, origin=[0,0], angle=0)
    return Ellipse(ep, mu, radius1, radius2, ShapeTransform(origin,angle))
end


"""
TODO
"""
function is_inside(ellipse::Ellipse, x::Real, y::Real)
    x, y = ellipse.transform \ [x; y]
    return (x/ellipse.radius1)^2 + (y/ellipse.radius2)^2 <= 1
end

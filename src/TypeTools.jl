"""
Functions to facilitate working with types.
"""
module TypeTools
export basetype
export parse_typedecl, parse_typesig, extract_typename, map_symbols
# export parameters, is_unionall_typevar, TypeMap
# export fields, properties, Property

#using MiscUtils
using MacroTools

# import Base.∘
# import Base.==
# import Base.fieldnames, Base.propertynames,  Base.getproperty
# import Base.show


"""
`is_unionall_typevar(tv, typ::Type)` returns true iff `tv` is a type variable
bound by a UnionAll in `typ`.
"""
function is_unionall_typevar(tv::TypeVar, typ::Type)
	while isa(typ, UnionAll)
		if tv == typ.var
			return true
		end
		typ = typ.body
	end
	false
end

is_unionall_typevar(::Any, ::Type) = false


"""
```
basetype(t::Type)
```
Returns a type withouts its parameters specified.
"""
# TODO: This doesn't give the obivous result for NTuple -- it returns Tuple.
# Desirable or not?
basetype(dt::DataType) = dt.name.wrapper
basetype(ut::UnionAll) = basetype(ut.body)

# Is this really needed?
# """
# 	freelastparameter(::Type)
#
# Free the last parameter of a type.  If it is already free, return the type unchanged.
#
# # Example
# ```
# julia> freelastparameter(Array{Int, 3})
# Array{Int, N} where N
# ```
# """
# function freelastparameter(::Type{A}) where A
# 	if isa(A, UnionAll)
# 		return A
# 	else
# 		basetype = A.name.wrapper
# 		return basetype{A.parameters[1:end-1]...}
# 	end
# end


# Parse type declarations (expressions of the form A{...} or A{...} <: B{...})
function parse_typedecl(head)
	@capture(head, (typesig_ <: supersig_) | typesig_) ||
	error("Header must be of the form `type` or `type <: type`")
	if supersig === nothing
		supersig = :( Any )
	end
	return (typesig, supersig)
end



# Parses a type signature, returning (name, parameters)
"""
	parse_typesig(expr) = (name::Symbol, params::[], qualparams::[])

Parse a type signature, e.g. an expression such as `A{T<:Number, Int}`. `name` is the
name of the type, `params` are the type parameters without qualifications, and `qualparams`
are the type parameters with any qualifications they may have.

# Example
```
parse_typesig(:(A{T<:Number,2})) = (:A, [:T, 2], [:(T<Number), 2])
```
"""
function parse_typesig(expr)
	qualparams = []
	if isa(expr, Symbol)
		typename = expr
	elseif expr.head == :curly
		typename = expr.args[1]
		append!(qualparams, expr.args[2:end])
	else
		error("$expr is not a valid type signature")
	end
	# Strip any qualifications that may be present on the parameters
	params = collect(map(qualparams) do qp
		if @capture(qp, param_ <: bound_)
			return param
		elseif @capture(qp, <: bound_)
			return gensym()
		else
			return qp
			end
	end)
	# for (i,qp) in enumerate(qualparams)
	# 	if @capture(qp, param_ <: bound_)
	# 		push!(params, param)
	# 	elseif @capture(qp, <: bound_)
	# 		push!(params, gensym())
	# 	else
	# 		push!(params, qp)
	# 	end
	# end
	(typename, params, qualparams)
end


function extract_typename(typesig)
    if isa(typesig, Symbol)
        typename = typesig
    elseif typesig.head == :curly
        typename = typesig.args[1]
    else
        error("Expected a type signature, got $typsig")
    end
end

# # This does not really extract type variables, but parameters without bounds
# function extract_typevars(typesig)
#     (tname, tparams) = parse_typeexpr(typesig)
#     typevars = Any[]
#     for (i,p) in enumerate(tparams)
#         if @capture(p, typevar_ <: bound_)
#             push!(typevars, typevar)
#         elseif @capture(p, <: bound_)
#             push!(typevars, gensym())
#        # elseif p isa Symbol
# 	 else
# 			push!(typevars, p)
# 		end
#     end
#     return typevars
# end


# Map select symbols within an expression.
# Wherever symbol `from_syms[i]` appears in `expr`, replace it with `to_exprs[i]`.
# Examples:
# map_symbols( :( x::Vector{T} ), [:T], [Float64]) = :( x::Vector{Float64} )
#
# SubType{S,T} <: SuperType{'x', T, Int64}
# map_smyobls( :( SuperType{'x', T, Int64} ), {S,T}, {A,B}) = :( SuperType{'x', B, Int64} )
function map_symbols(expr::Union{Expr, Symbol}, from_syms::Vector{Symbol}, to_exprs::Vector)
    new_expr = MacroTools.prewalk(expr) do s
        # print("Matching expression $s ...")
        if isa(s,Symbol)
            i = findfirst(s .== from_syms)
            if i == nothing
                # println("no match.")
                return s
            else
                # println("Mapping to $(to_syms[i])")
                return to_exprs[i]
            end
        else
            # println("no match")
            return s
        end
    end
    return new_expr
end

map_symbols(::Nothing, from_syms, to_syms) = nothing
# Anything that isn't an expression just gets passed through
map_symbols(x, from_syms, to_syms) = x


#=
The following stuff facilitates the propagation of type parameters to fields and supertypes.
=#


# Define fallbacks to make fieldnames and propertynames more consistent
#fieldnames(x) = fieldnames(typeof(x))
#propertynames(t::DataType) = fieldnames(t)
#propertynames(t::UnionAll) = fieldnames(t)
#propertynames(t::Type{<:Tuple}) = fieldnames(t)

# struct Property{T}
# 	name::Symbol
# 	typ::Type{T}
# 	default::Optional{T}
# end
# Property(n, t) = Property(n, t, nothing)
#
# propertytypes(x) = fieldtypes(typeof(x))
# propertytypes(x, private) = propertytypes(x)
#
# fields(x) = [Property(a,b) for (a,b) in zip(fieldnames(x), fieldtypes(x))]
# properties(x) = [Property(a,b) for (a,b) in zip(propertynames(x), propertytypes(x))]
#
# function show(io::IO, p::Property)
# 	print(io, p.name, "::")
# 	show(io, p.typ)
# 	if p.default != nothing
# 		print(io, " = ")
# 		show(io, p.default)
# 	end
# end
#


# NOT NEEDED.  Already provided by Base.nameof
#"""
#`name(::Type)` returns the name of a declared type.
#"""
#name(t::DataType) = t.name.name
#name(t::UnionAll) = name(t.body)


# """
# ```
# parameters(::Type)
# ```
# Returns the parameters of a declared type as a vector of types and/or typevars.
# """
# parameters(t::DataType) = t.parameters
# parameters(t::UnionAll) = parameters(t.body)

#
# """
# A `TypeMap` relates the formal parameters of a declared type A to the parameters
# of a related type B. (For example, B is a supertype of A or the type of a field
# of A.) The TypeMap can then be used to translate A{...} to the corresponding
# parameterized type B{...}.
# """
# struct TypeMap
# 	type_in::Type
# 	type_out::Type
# 	params::Vector
# end
#
#
# #  A Parameter is used in a TypeMap to represent a parameter to be looked up in the source.
# struct Param
# 	index::Int
# end

# This doesn't really have a plausible use case.
# """
# `TypeMap(T1, T2)` constructs a TypeMap from two formal types.
# """
# # Constructs a TypeMap from two formal types
# function TypeMap(type1::Type, type2::Type)
# 	# Parse the type signatures
# 	# Map parameters in type2 to those in type1
# 	params1 = parameters(type1)
# 	params2 = parameters(type2)
#
# 	pdict = Dict(zip(params1, 1:length(params1)))
# 	params = Vector()
# 	for p in params2
# 		if haskey(pdict, p)
# 			# Treat it as a reference to a member of params1
# 			push!(params, Param(pdict[p]))
# 		else
# 			# Treat it as a literal
# 			push!(params, p)
# 		end
# 	end
# 	TypeMap(basetype(type1), basetype(type2), params)
# end

# """
# `TypeMap(expr1, expr2)` creates a TypeMap from two expression, the first being
# the signature of a declared type and the second expressing a related type, whose
# parameters generally depend on those of the first. All types appearing in these
# expressions must be in scope.
# ```
# julia> tm = TypeMap(:(A{T,N}), :(B{N, Real}))
# julia> tm(A{Int, 3})
# B{3, Real}
# ```
# """
# # Constructs a TypeMap from two type signatues
# function TypeMap(expr1, expr2)
# 	# Parse the type signatures
# 	(tname1, params1) = parse_typeexpr(expr1)
# 	type1 = eval(current_module(), tname1)
#
# 	(tname2, params2) = parse_typeexpr(expr2)
# 	type2 = eval(current_module(), tname2)
#
# 	# Note, type1 and type2 will be as declared, since we eval'd just their names.
#
# 	# Map parameters in expr2 to those in expr1
# 	pdict = Dict(zip(params1, 1:length(params1)))
# 	params = Vector()
# 	for p in params2
# 		if haskey(pdict, p)
# 			# Treat it as a reference to a member of params1
# 			push!(params, Param(pdict[p]))
# 		else
# 			# Treat it as a literal
# 			push!(params, eval(current_module(), p))
# 		end
# 	end
# 	TypeMap(type1, type2, params)
# end
#
#
# """
# TypeMap(t::Type) creates the identity TypeMap for type T.
# """
# function TypeMap(t::Type)
# 	t_ad = basetype(t)
# 	nparams = length(parameters(t_ad));
# 	params = [Param(i) for i = 1:nparams]
# 	TypeMap(t_ad, t_ad, params)
# end
#
#
# """
# `(tm::TypeMap)(T::Type)` translates a parameterized type `T` as specified by
# type map `tm`. An error results if `T` is not of the right type for `tm`.
# """
# # TODO: Should this be a separate function?
# function (tmap::TypeMap)(typ::Type)
# 	@assert basetype(typ) == tmap.type_in "Expected input type $(tmap.type_in), got $typ"
#
# 	if length(tmap.params) == 0
# 		tmap.type_out
# 	else
# 		# collect bound TypeVars
# 		t_ = typ
# 		bound_vars = Set(TypeVar[])
# 		while isa(t_, UnionAll)
# 			push!(bound_vars, t_.var)
# 			t_ = t_.body
# 		end
#
# 		# map the parameters
# 		in_params = parameters(typ)
# 		out_params = Vector()
# 		# Loop over the TypeMap's params
# 		for a in tmap.params
# 			if isa(a, Param)
# 				# The parameter is formal (determined by the source type).
# 				push!(out_params, in_params[a.index])
# 			else
# 				# The parameter is actual.
# 				push!(out_params, a)
# 			end
# 		end
#
# 		type_out = tmap.type_out{out_params...}
#
# 		# Wrap bound TypeVars in type_out in UnionAll, prserving the order of
# 		# nesting that was in typ
# 		for i = length(out_params):-1:1
# 			if out_params[i] in bound_vars
# 				type_out = UnionAll(out_params[i], type_out)
# 			end
# 		end
# 		type_out
# 	end
# end
#
#
# ==(tm1::TypeMap, tm2::TypeMap) =
# 	tm1.type_in == tm2.type_in &&
# 	tm1.type_out == tm2.type_out &&
# 	tm1.params == tm2.params
#
#
# # Composes two TypeMaps
# function ∘(tmap2::TypeMap, tmap1::TypeMap)
# 	@assert tmap1.type_out == tmap2.type_in "Output type $(tmap1.type_out) of first TypeMap does not match input type $(tmap2.type_out) of second TypeMap"
#
# 	params1 = tmap1.params;
# 	params2 = Vector{Any}()
# 	for a in tmap2.params
# 		if isa(a, Param)
# 			push!(params2, params1[a.index])
# 		else
# 			push!(params2, a)
# 		end
# 	end
# 	TypeMap(tmap1.type_in, tmap2.type_out, params2)
# end

end

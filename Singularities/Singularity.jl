using  Singular, Polymake, LinearAlgebra
import  Nemo, GAP, Hecke
include("VecField.jl")

using Markdown

###############################################################################
#
#   Singularity Theory related functions TODO: Push into Singular.jl
#
###############################################################################

@doc Markdown.doc"""
    isisolated_singularity(f::spoly)
> Checks, whether the polynomial $f$ defines an isolated singularity.
> The ordering has to be local.
"""

function isisolated_singularity(f::Singular.spoly)
   R = parent(f)
   if !has_local_ordering(R)
      error("The ordering is not local.")
   end
   I1 = std(Singular.jacobian_ideal(f))
   if !iszerodim(I1)
      return false
   else
      return true
   end
end

@doc Markdown.doc"""
 milnor(f::Union{Singular.spoly{T}, Singular.sideal{T}})where T <: AbstractAlgebra.RingElem
> Computes the Milnor number, in case $f$ defines isolated hypersurface
> singularity or in case $f$ is the jacobian ideal of an isolated hypersurface singularity.
> The ordering has to be local.
"""

function milnor(f::Union{Singular.spoly{T}, Singular.sideal{T}})where T <: AbstractAlgebra.RingElem
   if typeof(f)==Singular.spoly{T}
    R = Singular.parent(f)
   elseif typeof(f)==Singular.sideal{T}
    R = Singular.base_ring(f)
   end

   if !has_local_ordering(R)
      error("The ordering is not local.")
   end
   if typeof(f)==Singular.spoly{T}
       J = std(Singular.jacobian_ideal(f))
   elseif typeof(f)==Singular.sideal{T}
       J = std(f)
   end

   if !iszerodim(J)
       error("The object does not define an isolated singularity.")
   else
       return vdim(std(J))
   end
end

@doc Markdown.doc"""
 tjurina(f::Union{Singular.spoly{T}, Singular.sideal{T}})where T <: AbstractAlgebra.RingElem
> Computes the Tjurina number, in case $f$ defines isolated hypersurface
> singularity or in case $f$ is the Tjurina ideal of an isolated hypersurface singularity.
> The ordering has to be local.
"""

function tjurina(f::Union{Singular.spoly{T}, Singular.sideal{T}})where T <: AbstractAlgebra.RingElem
   R = parent(f)
   if !has_local_ordering(R)
      error("The ordering is not local.")
   end
   I1 = std(Singular.jacobi(f))
   if !iszerodim(I1)
      error("The polynomial does not define an isolated singularity.")
   else
      I2 = Ideal(R, f)
      J = I1 + I2
    return vdim(std(J))
   end
end


@doc Markdown.doc"""
 kdeterminacy(f::spoly)
> Computes the Tjurina number, in case $f$ defines isolated hypersurface
> singularity. The ordering has to be local.
"""

function kdeterminacy(f::Singular.spoly)
   R = parent(f)
   if !has_local_ordering(R)
      error("The ordering is not local.")
   end
   I0 = std(Singular.jacobi(f))
   if !iszerodim(I0)
      error("The polynomial does not define an isolated singularity.")
   else

      I1 = std(MaximalIdeal(R, 1) * I0)
      I2 = std(MaximalIdeal(R, 2) * I0)

      n0 = total_degree(highcorner(I0)) + 2
      n1 = total_degree(highcorner(I1)) + 1
      n2 = total_degree(highcorner(I2))

      return min(n0, n1, n2)
   end
end

###############################################################################
#
#   Reconstruction algorithm
#
###############################################################################

function reconstruction_internal(I::Singular.sideal{Singular.spoly{Singular.n_Q}}, debug::Bool = false)
   # Compute GB if needed
   if !I.isGB
      I = std(I)
   end

   # Check if I is zero-dimensional
   !Singular.iszerodim(I) && error("Ideal is not zero-dimensional.")

   # Check if the ring has a local ordering
   R = base_ring(I)
   !has_local_ordering(R) && error("Ring needs local ordering.")

   # Compute Invariants
   debug && print("Compute Invariants:", "\n")
   n = nvars(R)
   mu = Int(milnor(I))
   bnd = Int(ceil(wdeg_bound_factor(n) * mu))
   debug && print("Number of Variables:", n, "\n")
   debug && print("Milnor number:", mu, "\n")
   debug && print("Weighted degree bound:", bnd, "\n")
   debug && print("DONE", "\n")

   # Compute logarithmic derivations
   debug && print("Compute logarithmic derivations:", "\n")
   Imin = Ideal(R, minimal_generating_set(I))
   Der = derivation_module(Imin)
   debug && print("Module of Logarithmic Derivations:", "\n", Der, "\n")
   debug && print("DONE", "\n")

   # Compute toral Lie algebra and coordinate change
   debug && print("Compute the new base ring and coordinate change:", "\n")
   B = get_sorted_kbasis(I)
   bnd2 = maximum(Singular.total_degree, B)
   #print("\n", "B=:", B,"\n")
   T = derivation_toral_lie_algebra_representation(Der, B)
   (S, phi, A) = get_coord_change(T, n, B)
   ##A = get_weight_vectors(Der)
   debug && print("Ring:", S, "\n", "Map:", "\n",phi, "\n")
   debug && print("DONE", "\n")

   # Change base ring of I
   debug && print("Computations in the new ring: ", "\n")
   Inew = Singular.Ideal(S, [change_base_ring(base_ring(S), I[i]; parent = S) for i in 1:Singular.ngens(I)])

   # Transform Inew into new coordinates
   phiinv = invertAlgebraMorphism(phi, bnd2)
   debug && print("phiinv=", "\n", phiinv, "\n")
   Inew = std(phiinv(Inew))
   #print("dernew=", "\n", derivation_module(Inew), "\n")
   soc = highcorner(Inew)
   exp = Singular.lead_exponent(soc)
   debug && print("Ideal:", Inew, "\n", "Socle: ", "\n", soc, "\n")
   debug && print("DONE", "\n")

   # Polytope computations to obtain candidate weights
   debug && print("Get the polytopes:", "\n")
   P = get_weight_polytope(A, bnd, exp)
   L = get_reduced_lattice_points(P, mu, exp)
   debug && print("The Polytope is given by the equations: ", P.EQUATIONS, "\n", "and inequalities:", P.INEQUALITIES, "\n", "reduced lattice points: ", "\n", L, "\n")
   debug && print("DONE", "\n")

   ## Candidate try outs:
   debug && print("Let the tryouts begin ...", "\n")
   s = size(L)
   # Get a minimal multihomogeneous generating system for Inew
   Imin = Singular.Ideal(S, minimal_weighted_homogeneous_generating_set(Inew, L[:, 2:s[2]]))
   for i in 1:s[1]
      Mon = get_weighted_homogeneous_monomials(S, L[i, :])
      Sp = build_candidate_space(S, Mon, Imin)
      for j in 1:10
         p = build_candidate_poly(S, Mon, Sp)
         if isisolated_singularity(p)
            debug && print("Result: ", phi(p), "\n")
            debug && print("Morphism: ", phi, "\n")
            return (p, phi)
         end
      end
   end

return false
end

###############################################################################
#
#   Auxilliary functions for Singular
#
###############################################################################

###############################################################################
#   Singular module from Array
###############################################################################

function Module(A::Array{Singular.spoly{T}, 2}
        ) where T <: AbstractAlgebra.RingElem
    R = Singular.parent(A[1])
    all(a -> Singular.parent(a) === R, A) || error("incompatible rings")
    cols = [ Singular.vector(R, A[:, i]...) for i in 1:size(A, 2) ]
    return Singular.Module(R, cols...)
end

###############################################################################
#   Concatenation of matrices
###############################################################################

@doc Markdown.doc"""
   hcat(A::smatrix{T}, B::smatrix{T}) where T <: AbstractAlgebra.RingElem
> Return the horizontal concatenation of $A$ and $B$.
> Assumes that the number of rows is the same in $A$ and $B$.
"""
function hcat(A::Singular.smatrix{T}, B::Singular.smatrix{T}) where T <: AbstractAlgebra.RingElem
   nr = Singular.nrows(A)
   R = base_ring(A)
   (Singular.nrows(B) != nr) && error("Matrices must have same number of rows.")
   (R != base_ring(B)) && error("Matrices are not over the same ring")

   nca = Singular.ncols(A)
   ncb = Singular.ncols(B)

   Z = Singular.zero_matrix(R, nr, nca + ncb)
   for i in 1:nr
      for j in 1:nca
         Z[i, j] = A[i, j]
      end

      for j in 1:ncb
         Z[i, nca + j] = B[i, j]
      end
   end
   return Z
end

@doc Markdown.doc"""
   hcat(A::Vector{ <: smatrix{T}) where T <: AbstractAlgebra.RingElem
> Return the horizontal concatenation of the matrices $A$.
> All component matrices need to have the same base ring and number of rows.
"""
function hcat(A::Vector{ <: Singular.smatrix{T} }) where T <: AbstractAlgebra.RingElem
  R = base_ring(A[1])
  nr = Singular.nrows(A[1])
  !all( M -> base_ring(M) == R, A) && error("Matrices are not over the same ring")
  !all( M -> Singular.nrows(M) == nr, A) && error("Matrices must have same number of rows.")

  l = length(A)
  nc = sum(Singular.ncols.(A))
  m = 0
  Z = Singular.zero_matrix(R, nr, nc)
  for i in 1:l
     nca = Singular.ncols(A[i])
     for j in m + 1:m + nca
        for k in 1:nr
           Z[k, j] = A[i][k, j - m]
        end
     end
     m += nca
  end
  return Z
end

@doc Markdown.doc"""
   vcat(A::smatrix{T}, B::smatrix{T}) where T <: AbstractAlgebra.RingElem
> Return the vertical concatenation of $A$ and $B$.
> Assumes that the number of columns is the same in $A$ and $B$.
"""
function vcat(A::Singular.smatrix{T}, B::Singular.smatrix{T}) where T <: AbstractAlgebra.RingElem
   nc = Singular.ncols(A)
   R = base_ring(A)
   (Singular.ncols(B) != nc) && error("Matrices must have same number of columns.")
   (R != base_ring(B)) && error("Matrices are not over the same ring")

   nra = Singular.nrows(A)
   nrb = Singular.nrows(B)

   Z = Singular.zero_matrix(R, nra + nrb, nc)
   for i in 1:nc
      for j in 1:nra
         Z[j ,i] = A[j, i]
      end

      for j in 1:nrb
         Z[nra + j, i] = B[j, i]
      end
   end
   return Z
end

@doc Markdown.doc"""
   vcat(A::Vector{ <: smatrix{T} })
> Return the vertical concatenation of the matrices $A$.
> All component matrices need to have the same base ring and number of columns.
"""
function vcat(A::Vector{ <: Singular.smatrix{T} }) where T <: AbstractAlgebra.RingElem
  R = base_ring(A[1])
  nc = Singular.ncols(A[1])
  !all( M -> base_ring(M) == R, A) && error("Matrices are not over the same ring")
  !all( M -> Singular.ncols(M) == nc, A) && error("Matrices must have same number of columns.")

  l = length(A)
  nr = sum(Singular.nrows.(A))
  m = 0
  Z = Singular.zero_matrix(R, nr, nc)
  for i in 1:l
     nra = Singular.nrows(A[i])
     for j in m + 1:m + nra
        for k in 1:nc
           Z[j, k] = A[i][j - m, k]
        end
     end
     m += nra
  end
  return Z
end

###############################################################################
#
#   Auxilliary functions for the reconstrution
#
###############################################################################

###############################################################################
#   GAP related functions
###############################################################################

@doc Markdown.doc"""
 convert_singularQQ_to_juliaQQ(L::Array{Array{n_Q,2},1})
> Converts a list of matrices with SINGULAR rational numbers into a
list of matrices containing JULIA rational numbers.

"""
function convert_singularQQ_to_juliaQQ(L::Array{Array{Singular.n_Q,2},1})
	r = size(L[1])[1]
	c = size(L[1])[2]
	LL = [[[Rational{BigInt}.(L[i][j,k]) for k in 1:c] for j in 1:r ] for i in 1:length(L)]
	return LL
end

##The following functions need rational SINGULAR coefficients
@doc Markdown.doc"""
 get_weight_lie_algebra(f::Union{Singular.spoly{Singular.n_Q}, Singular.sideal{Singular.spoly{Singular.n_Q}}, Singular.smodule{Singular.spoly{n_Q}}})
> Returns a list of matrices over the HECKE rationals, such that the eigenvalues of each matrix define the weight vectors of the ideal (generated by) f.
"""
function get_weight_lie_algebra(f::Union{Singular.spoly{Singular.n_Q}, Singular.sideal{Singular.spoly{Singular.n_Q}}, Singular.smodule{Singular.spoly{Singular.n_Q}}})

   # Avoid unnecessary syzygie computation
   if typeof(f) != Singular.smodule{Singular.spoly{Singular.n_Q}}
      L = derivation_module_get_linear_part(derivation_module(f))
   else
      L = derivation_module_get_linear_part(f)
   end

	L = convert_singularQQ_to_juliaQQ(L) ## Convert matrices from SingularQQ to JuliaQQ

	## Pass Objects to GAP to perform Lie algebra computations

	LG = GAP.julia_to_gap(L, Val(true)) ## Val(true) needed for recursion

	LG = GAP.Globals.LieAlgebra(GAP.Globals.Rationals, LG)

	C = GAP.Globals.LieCenter(GAP.Globals.CartanSubalgebra(LG))

	## List of Basis vectors of C
	C = GAP.Globals.BasisVectors(GAP.Globals.Basis(C))

	## Get the semi-simple parts and get a basis for the toral Lie algebra
	C = [GAP.Globals.JordanDecomposition(GAP.Globals.List(C[i]))[1] for i in 1:length(C)]

	C = GAP.Globals.LieAlgebra(GAP.Globals.Rationals, GAP.julia_to_gap(C))

	C = GAP.Globals.Basis(C)

	## Convert GAP output into list of HECKE matrices
	r = length(C)
	s = GAP.Globals.NrRows(C[1])
	L = [Hecke.matrix(reshape([Hecke.QQ(GAP.gap_to_julia(C[i][j, k])) for j in 1:s for k in 1:s], s, s)) for i in 1:r]

	return L
end

@doc Markdown.doc"""
 isAnalyticallyMonomial(f::Union{Singular.spoly{Singular.n_Q}, Singular.sideal{spoly{Singular.n_Q}}})
> Returns true if the ideal (generated by) f is a monomial ideal, false else.
> The ordering has to be local.
"""
function isAnalyticallyMonomial(f::Union{Singular.spoly{Singular.n_Q}, Singular.sideal{Singular.spoly{Singular.n_Q}}})
    if typeof(f)==Singular.spoly{T}
        R = Singular.parent(f)
    elseif typeof(f)==Singular.sideal{T}
        R = Singular.base_ring(f)
    end

    if !has_local_ordering(R)
      error("The ordering is not local.")
    else
        L = get_weight_lie_algebra(f)
        n = length(L)

        if typeof(f) == Singular.spoly{Singular.n_Q}
            if n == nvars(Singular.parent(f))
                return true
            else
                return false
            end
        elseif typeof(f) == Singular.sideal{Singular.spoly{Singular.n_Q}}
            if n == nvars(Singular.base_ring(f))
                return true
            else
                return false
            end
        end
    end
end

@doc Markdown.doc"""
 get_weight_vectors(I::Union{Singular.spoly{Singular.n_Q}, Singular.sideal{Singular.spoly{Singular.n_Q}}})
> Returns all weight vectors of I.
"""
function  get_weight_vectors(f::Union{Singular.spoly{Singular.n_Q}, Singular.sideal{Singular.spoly{Singular.n_Q}}, Singular.smodule{Singular.spoly{Singular.n_Q}}})

   if typeof(f) == Singular.spoly{Singular.n_Q}
      R = parent(f)
   end
   R = base_ring(f)
	L = get_weight_lie_algebra(f)
   #TODO: simultaneously diagonalize + base field change + coordchange
   # pass to Hecke....
   #TODO: Reduced row echelon form
	r = Hecke.ncols(L[1])
   s = length(L)

   # Return matrix with rational entries
	L = [reshape(Rational.(L[i]), r, r) for i in 1:s]
   A = zeros(Rational, s, r)
   for i in 1:s
      for j in 1:r
         A[i,j] = L[i][j,j]
      end
   end
	return A
end

###############################################################################
#   Hecke related functions
###############################################################################

@doc Markdown.doc"""
 getSplittingField(L::Array{Nemo.fmpq_mat,1})
> Returns a numberfield K over which the the input matrices are triagonalizable.
"""
function getSplittingField(L::Array{Nemo.fmpq_mat,1})

	Qx = Hecke.PolynomialRing(Hecke.base_ring(L[1])) ## Need this field for the next step
	SL = [Hecke.charpoly(Qx[1],L[i]) for i in 1:length(L)]

	return Hecke.splitting_field(SL)
end

function wdeg_bound_factor(n::Int)
   pr = 1
   b = 1

   for i in 1:n
      pr = Hecke.next_prime(pr)
      b = b * pr // (pr -1)
   end
   return b
end

###############################################################################
#   Polymake related functions
###############################################################################

function get_weight_polytope(A::Array{Rational{BigInt}, 2}, degbound::Int64, exp::Array{Int64, 1})
   n = size(A)[2]

   ##Building Cones
   E = Matrix{Rational{BigInt}}(I, n, n)
   #print("\n", "Matrices: ", E, A, "\n")
   c1 = Polymake.@pm polytope.Cone(RAYS = A)
   c2 = Polymake.@pm polytope.Cone(RAYS = E)
   c = Polymake.@pm polytope.intersection(c1, c2)
   #print("\n", "Cone: ", c1.RAYS, "\n", c2.RAYS, "\n", c.RAYS, "\n")
   ## Building weight polytope

   # Inequality via degree of socle
   B = zeros(Rational, 1, n+1)
   B[2:n+1] = -(exp.+2)
   B[1] = n * degbound

   # Inequalities via cone
   M = c.INEQUALITIES
   s = size(M)

   C = zeros(Rational, s[1], s[2]+1)

   for i in 1:s[1]
      for j in 1:s[2]
         C[i, j+1] = M[i, j]
      end
   end

   #Equations via cone
   M = c.EQUATIONS
   s = size(M)

   D = zeros(Rational, s[1], s[2]+1)
   for i in 1:s[1]
      for j in 1:s[2]
         D[i, j+1] = M[i, j]
      end
   end

   ## Return the weight polytope
   return Polymake.@pm polytope.Polytope(INEQUALITIES = Base.vcat(B , C), EQUATIONS = D)
end

function get_reduced_lattice_points(P::Polymake.BigObjectAllocated, mu::Int64, exp::Array{Int64, 1})

   n = length(exp)
   B = zeros(BigInt, 1, n+1)

   #Compute lattice points
   L = Polymake.polytope.lattice_points(P)
   L = BigInt.(L)
   s = size(L)
   #print("\n", "first L", L, "\n")
   for i in 1:s[1]
      socdeg = compute_socle_degree(exp, L[i, 2:n+1])

      if socdeg == false
         L[i, 1] = 0
      else
         L[i, 1] = socdeg
      end
   end
   #print("\n", "second L", L, "\n")
   for i in 1:s[1]
      if !any(x -> x == 0, L[i, :]) && gcd(L[i, :]) == 1 && milnor_check(mu, L[i,1], L[i, 2:n+1])
         B = Base.vcat(B, transpose(L[i, :]))
      end
   end
   #print("\n", "B = ", B, "\n")
   return Int.(B[2:end, :])
end

function get_weighted_homogeneous_monomials(R::Singular.PolyRing, dw::Array{Int, 1})

   n = nvars(R)
   n + 1 != length(dw) && error("Array does not have the right length.")

   B = zeros(Rational, 1, n+1)
   C = zeros(Rational, n, 1)
   E = Matrix{Rational}(I, n, n)
   B[2:n+1] = dw[2:n+1]
   B[1] = -dw[1]

   #Build polytope
   p = Polymake.@pm polytope.Polytope(INEQUALITIES = Base.hcat(C , E), EQUATIONS = B)
   L = BigInt.(Polymake.polytope.lattice_points(p))
   L = Int.(L)

   #Build monomials
   s = size(L)
   M = typeof(gen(R,1))[]
   for i in 1:s[1]
      push!(M, build_monomial(R, L[i, 2:n + 1]))
   end
   return M
end

###############################################################################
#   Ideal related functions
###############################################################################

function weighted_homogeneous_components(f::Singular.spoly, A::Array{Int, 2})
   D1 = Dict{Array{Int, 1}, MPolyBuildCtx}()
   D2 = Dict{Array{Int, 1}, typeof(f)}()
   s = size(A)
   R = parent(f)
   p = Singular.deepcopy(f)

   ## case f constant
   if Singular.isconstant(f)
      return push!(D2, (zeros(Int, 1, s[2])[1, :] => f))
   end
   ## Catch wrong number of columns
   s[2] != nvars(R) && error("Number of columns do not match number of ring
 variables")
   ## Compute the decomposition

   while p.ptr.cpp_object != C_NULL
   # Compute multi-degree
      deg = zeros(Int, 1, s[1])[1, :]
      exp = Singular.lead_exponent(p)
      for i in 1:s[1]
         deg[i] = dot(exp, A[i, :])
      end

      # Add monomial to corresponding component
      if haskey(D1, deg)
         Singular.push_term!(D1[deg], lc(p), exp)
      else
         M = Singular.MPolyBuildCtx(R)
         Singular.push_term!(M, lc(p), exp)
         push!(D1, deg => M)
      end
      # Change head of p
      p.ptr = Singular.libSingular.pNext(p.ptr)
   end

   ## Return final dictionary
   for d in keys(D1)
      push!(D2, d => Singular.finish(D1[d]))
   end
   return D2
end

function minimal_weighted_homogeneous_generating_set(I::Singular.sideal,
           A::Array{Int, 2})
   R = base_ring(I)
   n = Singular.ngens(I)

   !has_local_ordering(R) && error("Ring needs local ordering.")

   L = Array{elem_type(R), 1}()
   for i in 1:n
      D = weighted_homogeneous_components(I[i], A)
      for k in keys(D)
         push!(L, D[k])
      end
   end

   return Singular.minimal_generating_set(Singular.Ideal(R, L))
end

###############################################################################
#   Morphism related functions
###############################################################################

function invertAlgebraMorphism(phi::Singular.SAlgHom, n::Int)
# assume domain = codomain
   Im = phi.image
   R = phi.codomain
   maxId = MaximalIdeal(R, n)
   n = nvars(R)
   S, = Singular.PolynomialRing(base_ring(R), ["x$i" for i in 1:2*n], ordering = :lex)

   # Build embedding and ideal
   emb = Singular.AlgebraHomomorphism(R, S, gens(MaximalIdeal(S, 1))[1:n])
   IdAr1 = [gen(S, n+ i) - emb(Im[i]) for i in 1:n]
   SiD = Singular.Ideal(S, IdAr1) + emb(maxId)

   # Compute inverse and project down
   SiD = Singular.std(SiD, complete_reduction = true)
   IdAr2 = [R(0) for i in 1:n]
   proj = Singular.AlgebraHomomorphism(S, R, Singular.vcat(IdAr2, gens(MaximalIdeal(R, 1))))

   return Singular.AlgebraHomomorphism(R, R, [proj(reduce(gen(S, i), SiD)) for i in 1:n])
end

###############################################################################
#   Computation of invariants
###############################################################################

function compute_socle_degree(exp::Array{Int64, 1}, w::Array{Integer, 1})
   res = (dot(exp, w) + 2 * sum(w)) / length(exp)
   if denominator(Rational{BigInt}(res)) == 1
      return BigInt(res)
   else
      return false
   end
end

function milnor_check(mu::Int64, d::BigInt, w::Array{Integer, 1})

   x = d // w[1] - 1
   for i in 2:length(w)
      x = x*(d // w[i] - 1)
   end
   return mu == x
end

function build_monomial(R::Singular.PolyRing, exp::Array{Int, 1})
   M = Singular.MPolyBuildCtx(R)
   return Singular.push_term!(M, base_ring(R)(1), exp)
end

###############################################################################
#   General auxilliary functions
###############################################################################

function build_syz_matrix(R::Singular.PolyRing, M::Array{T, 1}, I::Singular.sideal) where T<: Singular.RingElem
   k = length(M)
   n = Singular.nvars(R)
   m = Singular.ngens(I)

   IdMat1 = Singular.zero_matrix(R, 1, k)
   IdMat2 = Singular.zero_matrix(R, 1, m)
   for i in 1:k
      IdMat1[1, i] = M[i]
   end

   for i in 1:m
      IdMat2[1, i] = I[i]
   end

   Jac = Singular.jacobian_matrix(M)
   Mat = vcat(IdMat1, Singular.transpose(Jac))

   B = zeros(R, n + 1, k + (n+1)*m)

   # AbstractAlgebra unterstützt  B[1:n + 1, 1:k] = Mat[:,:] nicht ...
   for i in 1:n + 1
      for j in 1:k
         B[i, j] = Mat[i, j]
      end
   end
   # Fill other values
   for i in 1:n+1
      for j in k + (i - 1)*m + 1:k + (i-1)*m + m
         B[i, j] = IdMat2[1, j-k-(i-1)*m]
      end
   end

   return B
end

function build_candidate_space(R::Singular.PolyRing, M::Array{T, 1}, I::Singular.sideal) where T<: Singular.RingElem
   k = length(M)
   B = build_syz_matrix(R, M, I)
   Mod = syz(Module(B))
   L = Array(Mod[1])[1:k]
   for i in 2:Singular.ngens(Mod)
      L = Base.hcat(L, Array(Mod[i])[1:k])
   end
   #Remove zeros in L
   Mod = Module(L)
   Singular.libSingular.idSkipZeroes(Mod.ptr)

   # Return generators of degree 0 part of module
   Mod = jet(Mod, 0)
   S = Array(Mod[1])[1:k]
   for i in 2:Singular.ngens(Mod)
      S = Base.hcat(S, Array(Mod[i])[1:k])
   end
   return Array(transpose(S))
end

function build_candidate_poly(R::Singular.PolyRing, M::Array{T, 1}, S::Array{T, 2}) where T<: Singular.RingElem
   s = size(S)
   P = [dot(S[i, :], M) for i in 1:s[1]]
   Ran = R.(rand(-100:100, length(P)))
   return dot(Ran, P)
end


import AbstractAlgebra, Hecke, Nemo

using Markdown

################################################################################
#
# Diagonalization Stuff by Raul
#
################################################################################

@doc Markdown.doc"""
 multispectrum(L::Array{S, 1}; splitting_field::Bool=false)
> Returns a dictionary where the keys are arrays of eigenvalues of the matrices
> in $L$ and where the values are the corresponding multiplicities. If
> "splitting_field" is set to "true", the multispectrum is computed over a
> splitting field of the elements of $L$.
"""
function multispectrum(L::Array{S, 1}; splitting_field::Bool=false) where S <:Hecke.MatElem{T} where T<: Hecke.FieldElem

  n = length(L)
  s = Hecke.nrows(L[1])

  # Pass to splitting field
  if splitting_field
    K = getSplittingField(L)
    L = [Hecke.matrix(reshape([K(L[i][j, k]) for j in 1:s for k in 1:s], s, s)) for i in 1:n]
  else
    K = parent(L[1][1, 1])
  end

  # Compute eigenvalues and their multiplicities
  CE = common_eigenspaces(list_of_eigenspaces(L))
  D = Dict{Array{typeof(L[1][1, 1]), 1}, Int64}()

  for (k, v) in CE
    push!(D, k => Hecke.nrows(v))
  end

  return D
end

@doc Markdown.doc"""
 multieigenspace(L::Array{S, 1}; splitting_field::Bool=false)
> Returns a a matrix, whose rows contain the eigenvectors the multi-eigenvalue.
> 
"""
function multieigenspace(L::Array{S, 1}, E::Array{T, 1}) where S <:Hecke.MatElem{T} where T<: Hecke.FieldElem

   n = length(L)
   s = Hecke.nrows(L[1])
   K = parent(L[1][1, 1])
   id = Hecke.identity_matrix(K, s)

   n != length(E) && error("Array dimension mismatch.")  
   
   # Compute eigenspaces
   Mat = L[1] - E[1]*id
   for i in 2:n
      Mat = Nemo.vcat(Mat, L[i] - E[i]*id)
   end
   return transpose(Hecke.kernel(Mat)[2])
end

@doc Markdown.doc"""
 simuldiag(L::Array{S, 1}; splitting_field::Bool=false)
> Returns a tuple whose first entry is the transformation matrix and whose
> second entry is an array of matrices containing the diagonal forms of
> the elements of $L$. If "splitting_field" is set to "true", the
> computation is performed over a splitting field of the elements of $L$.
> Default value is "false".
> If "check_commute_pairwise" is set to true, the algorithm checks, whether
> the matrices in $L$ commute pairwise. Default value is "true".
"""
function simuldiag(L::Array{S, 1}; splitting_field::Bool=false, check_commute_pairwise::Bool=true) where S <:Hecke.MatElem{T} where T<: Hecke.FieldElem

  if check_commute_pairwise
    if !commute_pairwise(L)
      error("Matrices do not commute pairwise.")
    end
  end

  n = length(L)
  s = Hecke.nrows(L[1])

  # Pass to splitting field
  if splitting_field
    K = getSplittingField(L)
    L = [Hecke.matrix(reshape([K(L[i][j, k]) for j in 1:s for k in 1:s], s, s)) for i in 1:n]
  else
    K = parent(L[1][1,1])
  end

  # Compute transformation marix
  CE = common_eigenspaces(list_of_eigenspaces(L))
  A =  Hecke.vcat(collect(values(CE)))

  # Compute diagonal forms
  D = [AbstractAlgebra.zero_matrix(K, s, s) for i in 1:n]
  m = 0
  for (v, k) in CE
    nr = Hecke.nrows(k)
    for j in 1:nr
      for i in 1:n
        D[i][m + j, m + j] = v[i]
      end
    end
    m = m + nr
  end
  return (A, D)
end

#################################################################################
#
# Possibly useful functions
#
#################################################################################

@doc Markdown.doc"""
 intersect_spaces(A::Hecke.MatElem{T}, B::Hecke.MatElem{T})
> Given two subvectorspaces $U$ and $V,$ whose bases are the rows of $A$ and $B$,
> this function returns a matrix whose rows are a basis of $U \cap V.$
"""
function intersect_spaces(A::Hecke.MatElem{T}, B::Hecke.MatElem{T}) where T <: Hecke.FieldElem

  n = Hecke.nrows(A)
  M = Hecke.vcat(A, B)
  N = Hecke.kernel(M, side=:left)[2]
  return N[:, 1:n]*A
end

@doc Markdown.doc"""
 getSplittingField(L::Array{S, 1})
> Given an array $L$ of matrices, this function returns a splitting field of the
> matrices in $L$.
"""
function getSplittingField(L::Array{S, 1}) where S <: Hecke.MatElem{T} where T <: Hecke.FieldElem

	Qx = Hecke.PolynomialRing(Hecke.base_ring(L[1])) # Need this field for the next step
	SL = [Hecke.charpoly(Qx[1],L[i]) for i in 1:length(L)]
	return Hecke.splitting_field(SL)
end


#################################################################################
#
# Auxilliary functions
#
#################################################################################

function commute_pairwise(L::Array{S, 1}) where S <: Hecke.MatElem{T} where T <: Hecke.FieldElem

  n = length(L)

  # Check pairwise commutativity
  for i in 1:n
    for j in i:n
      if L[i]*L[j]!=L[j]*L[i]
        return false
      end
    end
  end
  return true
end

function eigenspaces(M::Hecke.MatElem{T}) where T<: Hecke.FieldElem

  S = Hecke.spectrum(M)
  nr = Hecke.nrows(M)

  # Check if charpoly(M) factorizes into linear factors
  if sum(collect(values(S))) != nr
    error("Matrix is not diagonalizable.")
  end

  n = length(S)
  L = Dict{Array{typeof(M[1, 1]), 1}, typeof(M)}()
  for k in keys(S)
    push!(L, [k] => Hecke.vcat(Hecke.eigenspace(M, k)))
  end

  # Check if the eigenspaces of M form a basis of the vector space.
  if sum(Hecke.nrows.(collect(values(L)))) != nr
    error("Matrix is not diagonalizable.")
  end
  return L
end

function intersect_eigenspaces(L1::Dict{Array{T, 1}, S}, L2::Dict{Array{T, 1}, S}) where S<:Hecke.MatElem{T} where T <: Hecke.FieldElem

  L = Dict{keytype(L1), valtype(L1)}()
  for (k1, v1) in L1
    for (k2, v2) in L2
      I = intersect_spaces(v1, v2)
      if !iszero(I)
        push!(L, Nemo.vcat(k1, k2)  => I)
      end
    end
  end
  return L
end

function list_of_eigenspaces(L::Array{S, 1}) where S <: Hecke.MatElem{T} where T <: Hecke.FieldElem

  n = length(L)
  LL = [eigenspaces(L[1])]
  for i in 2:n
    push!(LL, eigenspaces(L[i]))
  end
  return LL
end

function common_eigenspaces(L::Array{Dict{Array{T, 1}, S}, 1}) where S<:Hecke.MatElem{T} where T<:Hecke.FieldElem

  n = length(L)
  if n==1
    return L[1]
  end
  k = BigInt(floor(n/2))
  return intersect_eigenspaces(common_eigenspaces(L[1:k]), common_eigenspaces(L[k+1:n]))
end


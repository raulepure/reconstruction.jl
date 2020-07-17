module Reconstruction
using Oscar
import Markdown, Singular, GAP, Hecke, Polymake

include("/Singularities/Diagonal.jl")
include("/Singularities/Singularity.jl")
include("/Singularities/VecField.jl")


export reconstruction

function reconstruction(I::Oscar.MPolyIdeal; debug::Bool = false)
   Oscar.singular_assure(I)
   B = base_ring(I)
   J = changeOrderOfBasering(I.gens.S, :negdegrevlex)
   return reconstruction_internal(J, debug)
end

#changes order of basering to lex via an algebra homomorphism
function changeOrderOfBasering(I::Singular.sideal, ordering::Symbol = lex)
    R = I.base_ring;
    G = Singular.gens(I.base_ring);
    Gstrich = string.(G);
    S, G = Singular.PolynomialRing(R.base_ring, Gstrich, ordering = ordering)
    res = Singular.AlgebraHomomorphism(R, S, G);
    return res(I)
end


end
using .Reconstruction

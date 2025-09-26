hankel1(n,x) = besselj(n,x) + 1im*bessely(n,x)
hankel1_deriv(n,x) = (hankel1(n-1,x) - hankel1(n+1,x))/2
besseli_deriv(n,x) = (besseli(n-1,x) + besseli(n+1,x))/2
besselj_deriv(n,x) = (besselj(n-1,x) - besselj(n+1,x))/2
besselk_deriv(n,x) = -(besselk(n-1,x) + besselk(n+1,x))/2

struct Fiber
    l::Int64
    n1::Float64 
    n2::Float64
    a::Float64
    k::Float64
    u::Float64
    v::Float64
    w::Float64
end

function fiber_cond(x,p)
    l,n1,n2,v = p    
    u,w = x
    
    f1 = (besselj_deriv(l,u)/(u*besselj(l,u)) + besselk_deriv(l,w)/(w*besselk(l,w)))
    f2 = (n1^2*besselj_deriv(l,u)/(u*besselj(l,u)) + n2^2*besselk_deriv(l,w)/(w*besselk(l,w)))
    f3 = (1/u^2 + 1/w^2)*(n1^2/u^2 + n2^2/w^2)
    
    cond1 = f1*f2 - l^2*f3
    cond2 = u^2 + w^2 - v^2
    
    return [cond1, cond2]
end

function fiber_parameters(; l::Int64 = 1, n1::Real = 1.4537, n2::Real = 1., a::Real = 200e-9, λ::Real = 852.3e-9)

    k = 2*π/λ 

    v = a*k*sqrt(n1^2 - n2^2)
    u0 = [2.,1.]
    prob = NonlinearProblem(fiber_cond,u0,[l, n1, n2, v])
    sol = solve(prob)
    u, w = sol.u 
    
    return Fiber(l,n1,n2,a,k,u,v,w)
end

function fiber_coupling(fiber::Fiber, intensity::Matrix{Float64}, dx::Float64; Lcav::Float64 = 0.33, dipole_factor::Float64 = 1.0)
    α = 1/137
    μ = 4.4837
    
    Vmode = Lcav*fiber.a^2*sum(intensity)*dx*dx
    
    # in units of 2π x MHz
    g0 = dipole_factor*μ*c*a0*sqrt(α/Vmode)*sqrt(fiber.k/(2*π))
    g0 /= 1e6

    return g0
end

function fiber_mode(pos::Vector{Float64}, fiber::Fiber; coordinates::Symbol = :cartesian, circ::Int64 = 1, t::Real = 0)

    l = fiber.l; n1 = fiber.n1; n2 = fiber.n2; a = fiber.a; 
    k = fiber.k; u = fiber.u; v = fiber.v; w = fiber.w;

    ω = k*c
    
    if coordinates == :cartesian
        x,y,z = pos
        r = norm(pos)
        ϕ = atan(y,x) 
    else
        r,ϕ,z = pos 
    end

    β = k*sqrt((n1^2/u^2 + n2^2/w^2)/(1/u^2 + 1/w^2))     
    s = l*(1/u^2 + 1/w^2)/(besselj_deriv(l,u)/(u*besselj(l,u)) + besselk_deriv(l,w)/(w*besselk(l,w)));
            
    if r >= 1
        er = 1im*(a*β/w)*(besselj(l,u)/besselk(l,w))*(besselk_deriv(l,w*r) - l*s*besselk(l,w*r)/(w*r))
        eϕ = -circ*(a*β/w)*(besselj(l,u)/besselk(l,w))*(l*besselk(l,w*r)/(w*r) - s*besselk_deriv(l,w*r))
        ez = besselk(l,w*r)*besselj(l,u)/besselk(l,w)
    else
        er = -1im*(a*β/u)*(besselj_deriv(l,u*r) - l*s*besselj(l,u*r)/(u*r))  
        eϕ = circ*(a*β/u)*(l*besselj(l,u*r)/(u*r) - s*besselj_deriv(l,u*r))
        ez = besselj(l,u*r)
    end
    
    field = [er,eϕ,ez]
    
    if coordinates == :cartesian
        ex = er*cos(ϕ) - eϕ*sin(ϕ)
        ey = er*sin(ϕ) + eϕ*cos(ϕ)
        field = [ex,ey,ez]    
    end
    
    field *= exp(1im*(ω*t + circ*l*ϕ - β*z))
    
    return field
end

function fiber_mode_linear(pos::Vector{Float64}, fiber::Fiber; pol_angle::Real = 0., coordinates::Symbol = :cartesian, circ::Int64 = 1) 

    E_left = fiber_mode(pos,fiber;coordinates=coordinates,circ=circ)
    E_right = fiber_mode(pos,fiber;coordinates=coordinates,circ=-circ)

    E_lin = (E_left*exp(-1im*deg2rad(pol_angle)) + E_right*exp(1im*deg2rad(pol_angle)))/sqrt(2)

    return E_lin
end
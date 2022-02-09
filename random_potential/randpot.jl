# Create 2d Random Potential

#module randpot
using FFTW
using Interpolations
using Statistics

#export Correlation, GaussCorr, RandPot2D,updatePot!

abstract type Correlation end

struct GaussCorr <: Correlation
    s::Float64 #sqrt(variance)
    l::Float64
    n::Array{Int32}
    c::Array{Float64}
    cf::Array{ComplexF64}
    fxf::Array{ComplexF64}
    fyf::Array{ComplexF64}
    ext::Array{Float64} #[x0,x1,y0,y1]
    kx::Array{Float64}
    ky::Array{Float64}
    L::Array{Float64}
    x::StepRangeLen
    y::StepRangeLen

    function GaussCorr(;s=1,l=0.05,n=[1024,1024],ext=[-0.5,0.5,-0.5,0.5])
        corr(x,y)=s*exp(-(x^2+y^2)/(l^2))
        L=zeros(Float64,2)
        L[1]=ext[2]-ext[1]
        L[2]=ext[4]-ext[3]
        x=range(ext[1],step=L[1]/n[1],length=n[1]+1)
        y=range(ext[3],step=L[2]/n[2],length=n[2]+1)
        ctmp=[corr(X,Y) for X in x[1:n[1]], Y in y[1:n[2]]]
        c=ifftshift(ctmp,(-Int(n[1]/2),-Int(n[2]/2)))
        kx=[ ((i>n[1]÷2) ? float(i-n[1]) : float(i))*2*π/L[1] for i in range(0,length=n[1]),j in range(0,length=n[2])]
        ky=[ ((j>n[2]÷2) ? float(j-n[2]) : float(j))*2*π/L[2] for i in range(0,length=n[1]), j in range(0,length=n[2])]
        cf=sqrt.(fft(c))
        fxf=1im*kx .* cf
        fyf=1im*ky .* cf
       return new(s,l,n,c,cf,fxf,fyf,ext,kx,ky,L,x,y)
   end
end

function gaussani(x,y,l;p=0)
    a=1/l[1]
    b=1/l[2]
    Mx=(a*cos(p))^2+(b*sin(p))^2
    My=(a*sin(p))^2+(b*cos(p))^2
    Mxy=2*(a^2-b^2)*sin(p)*cos(p)
    return exp(-(Mx*x*x+My*y*y+Mxy*x*y))
end

struct AnisoGaussCorr <: Correlation
    s::Float64 #sqrt(variance)
    l::Array{Float64}
    phi::Float64
    n::Array{Int32}
    c::Array{Float64}
    cf::Array{ComplexF64}
    fxf::Array{ComplexF64}
    fyf::Array{ComplexF64}
    ext::Array{Float64} #[x0,x1,y0,y1]
    kx::Array{Float64}
    ky::Array{Float64}
    L::Array{Float64}
    x::StepRangeLen
    y::StepRangeLen

    function AnisoGaussCorr(;s=1,l=[0.05,0.025],n=[1024,1024],ext=[-0.5,0.5,-0.5,0.5],phi=0)
        corr(x,y)=gaussani(x,y,l;p=phi)
        L=zeros(Float64,2)
        L[1]=ext[2]-ext[1]
        L[2]=ext[4]-ext[3]
        x=range(ext[1],step=L[1]/n[1],length=n[1]+1)
        y=range(ext[3],step=L[2]/n[2],length=n[2]+1)
        ctmp=[corr(X,Y) for X in x[1:n[1]], Y in y[1:n[2]]]
        c=ifftshift(ctmp,(-Int(n[1]/2),-Int(n[2]/2)))
        kx=[ ((i>n[1]÷2) ? float(i-n[1]) : float(i))*2*π/L[1] for i in range(0,length=n[1]),j in range(0,length=n[2])]
        ky=[ ((j>n[2]÷2) ? float(j-n[2]) : float(j))*2*π/L[2] for i in range(0,length=n[1]), j in range(0,length=n[2])]
        cf=sqrt.(fft(c))
        fxf=1im*kx .* cf
        fyf=1im*ky .* cf
       return new(s,l,phi,n,c,cf,fxf,fyf,ext,kx,ky,L,x,y)
   end
end


function randomphases!(phase)
    n=size(phase)
    for i in range(0,length=Int(n[1]/2))
        for j in range(0,length=Int(n[2]))
            tmp=rand()*2*π
            phase[i+1,j+1]=exp(1im*tmp)
            phase[n[1]-i,n[2]-j]=exp(-1im*tmp)
        end
    end
    phase[1,1]=0
    phase[1,Int(n[2]/2)]=1
    phase[Int(n[1]/2),1]=1
    phase[Int(n[1]/2),Int(n[2]/2)]=1
end

function udatePot2d!(c::Correlation,phases::Array{ComplexF64},u::Array{Float64},fx::Array{Float64},fy::Array{Float64})
    N=sqrt(2*c.n[1]*c.n[2])
    q=ifft(c.cf .* phases)
    u[1:c.n[1],1:c.n[2]]=N*real.(q)
    #println(maximum(imag(q)))
    fx[1:c.n[1],1:c.n[2]]=N*real.(ifft(c.fxf .* phases))
    fy[1:c.n[1],1:c.n[2]]=N*real.(ifft(c.fyf .* phases))
end

function createInterpolations(c::Correlation,u::Array{Float64},fx::Array{Float64},fy::Array{Float64})
    U = scale(interpolate(u, BSpline(Cubic(Line(OnGrid())))),c.x,c.y)
    Fx = scale(interpolate(fx, BSpline(Cubic(Line(OnGrid())))),c.x,c.y)
    Fy = scale(interpolate(fy, BSpline(Cubic(Line(OnGrid())))),c.x,c.y)
    return(U,Fx,Fy)
end

# function pbc{T<:Correlation}(x::Float64,c::T,i::Int64)
#     return (mod((x-c.ext[(i-1)*2+1]),c.L[i])+c.ext[(i-1)*2+1])
# end

mutable struct RandPot2D{T<:Correlation}
    corr::T
    phases::Array{ComplexF64}
    ui::Array{Float64}
    fxi::Array{Float64}
    fyi::Array{Float64}
    u::Array{Float64}
    fx::Array{Float64}
    fy::Array{Float64}
    Ui::ScaledInterpolation
    Fxi::ScaledInterpolation
    Fyi::ScaledInterpolation
    U::Function
    Fx::Function
    Fy::Function
    x::Array{Float64}
    y::Array{Float64}


    function RandPot2D{T}(c::T) where {T <: Correlation}
        phases=zeros(ComplexF64,c.n...)
        ui=zeros((c.n[1]+1,c.n[2]+1))
        fxi=zeros((c.n[1]+1,c.n[2]+1))
        fyi=zeros((c.n[1]+1,c.n[2]+1))

        randomphases!(phases)
        udatePot2d!(c,phases,ui,fxi,fyi)
        s=std(ui)
        ui=ui/s*c.s
        fxi=fxi/s*c.s
        fyi=fyi/s*c.s
        u=view(ui,1:c.n[1],1:c.n[2])
        fx=view(fxi,1:c.n[1],1:c.n[2])
        fy=view(fyi,1:c.n[1],1:c.n[2])
        x=collect(c.x[1:c.n[1]])
        y=collect(c.y[1:c.n[2]])
        x0=c.ext[1]
        y0=c.ext[3]
        lx=c.L[1]
        ly=c.L[2]
        Ui,Fxi,Fyi=createInterpolations(c,ui,fxi,fyi)
        U(x::Float64,y::Float64)=Ui(mod(x-x0,lx)+x0,mod(y-y0,ly)+y0)
        Fx(x::Float64,y::Float64)=Fxi(mod(x-x0,lx)+x0,mod(y-y0,ly)+y0)
        Fy(x::Float64,y::Float64)=Fyi(mod(x-x0,lx)+x0,mod(y-y0,ly)+y0)
        # U(x::Float64,y::Float64)=Ui(pbc(x,c,1),pbc(y,c,2))
        # Fx(x::Float64,y::Float64)=Fxi(pbc(x,c,1),pbc(y,c,2))
        # Fy(x::Float64,y::Float64)=Fyi(pbc(x,c,1),pbc(y,c,2))

        return(new{typeof(c)}(c,phases,ui,fxi,fyi,u,fx,fy,Ui,Fxi,Fyi,U,Fx,Fy,x,y))
    end
end

function updatePot!(pot::RandPot2D)
    c=pot.corr
    randomphases!(pot.phases)
    udatePot2d!(c,pot.phases,pot.ui,pot.fxi,pot.fyi)

    pot.u=view(pot.ui,1:c.n[1],1:c.n[2])
    pot.fx=view(pot.fxi,1:c.n[1],1:c.n[2])
    pot.fy=view(pot.fyi,1:c.n[1],1:c.n[2])
    s=std(pot.u)
    pot.ui=pot.ui/s
    pot.fxi=pot.fxi/s
    pot.fyi=pot.fyi/s

    Ui,Fxi,Fyi=createInterpolations(c,pot.ui,pot.fxi,pot.fyi)
    # U(x,y)=Ui(pbc(x,c,1),pbc(y,c,2))
    # Fx(x,y)=Fxi(pbc(x,c,1),pbc(y,c,2))
    # Fy(x,y)=Fyi(pbc(x,c,1),pbc(y,c,2))
    x0=c.ext[1]
    y0=c.ext[3]
    lx=c.L[1]
    ly=c.L[2]
    U(x::Float64,y::Float64)=Ui(mod(x-x0,lx)+x0,mod(y-y0,ly)+y0)
    Fx(x::Float64,y::Float64)=Fxi(mod(x-x0,lx)+x0,mod(y-y0,ly)+y0)
    Fy(x::Float64,y::Float64)=Fyi(mod(x-x0,lx)+x0,mod(y-y0,ly)+y0)
    pot.U=U
    pot.Fx=Fx
    pot.Fy=Fy
end #function updatePot!(args)

#end  # module randpot

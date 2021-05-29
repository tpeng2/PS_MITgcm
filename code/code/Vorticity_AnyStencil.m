function Omega=Vorticity_AnyStencil(dx,dy,U,V,npts)
    %Calculates the vorticity of a 2D (U,V) field in an UNIFORM GRID.
    %dx and dy are the grid spacings, U, V are the vector components
    %npts gives the stencil size (uses centered differences, so npts must
    %be odd. If npts is even, adds 1).
    
    %U and V grids in ***meshgrid*** format.
    
    %Fernando Zigunov, 2019
    if mod(npts,2)==0
        npts=npts+1;%If npts is even, adds 1 for centered difference
    end
    
    StencilX=FiniteDifferenceStencil(dx,npts,1);
    StencilY=FiniteDifferenceStencil(dy,npts,1).';
    
    %dudx=conv2(U,StencilX,'same');
    dudy=conv2(U,StencilY,'same'); %Calculates derivative with stencil by convolution
    
    dvdx=conv2(V,StencilX,'same');
    %dvdy=conv2(V,StencilY,'same');
    
    Omega=dudy-dvdx;    
end
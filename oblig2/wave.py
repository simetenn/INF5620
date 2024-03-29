import nose.tools as nt
import pylab as p
import sys


def solver(I, V, q, f ,b, Lx, dx, Ly, dy, T, dt, version="scalar"):
    """
    The solver function
    """
    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    u = p.zeros((Nx,Ny,Nt))

    x = p.linspace(0, Lx, Nx)
    y = p.linspace(0, Ly, Ny)

    #Create initial conditions stuff
    if version == 'vectorized':
        init_vec(u, I, V, x, y, dt)
        advance_vec(u, q, f, x, y, b, dx, dy, Nt, dt)
    else:
        init(u, I, V, Nx, Ny, x, y, dt)
        advance(u, q, f, x, y, b, Nx, dx, Ny, dy, Nt ,dt)
                
    return u


    
    
def init(u, I, V, Nx, Ny, x, y, dt):
    """
    Set initial conditions
    """
    #Initial position
    for i in xrange(0, Nx):
        for j in xrange(0, Ny):
            u[i,j,0] = I(x[i],y[j])
    #previous position
    for i in xrange(0, Nx):
        for j in xrange(0, Ny):
            u[i,j,-1] = u[i,j,0] - dt*V(x[i],y[j])

            

def advance(u, q, f, x, y, b, Nx, dx, Ny, dy, Nt, dt):
    """
    Advance in time
    """
    dt2 = dt**2
    dtdx2 = dt2/(2*dx**2)
    dtdy2 = dt2/(2*dy**2)
    
    for n in xrange(0,Nt-1):
        for i in xrange(1, Nx-1):
            for j in xrange(1, Ny-1):
                scheme(u, q, f,  i, j, n, i+1, i-1, j+1, j-1, x, y, dtdx2, dtdy2, dt2, dt, b)

            
        #Update boundaries
        for i in xrange(1,Nx-1):
            j = 0;
            scheme(u, q, f, i, j, n, i+1, i-1, j+1, j+1, x, y, dtdx2, dtdy2, dt2, dt, b)
            
            j = Ny-1;
            scheme(u, q, f, i, j, n, i+1, i-1, j-1, j-1, x, y, dtdx2, dtdy2, dt2, dt, b)
            
        for j in xrange(1,Ny-1):
            i = 0;
            scheme(u, q, f, i, j, n, i+1, i+1, j+1, j-1, x, y, dtdx2, dtdy2, dt2, dt, b)
            
            i = Nx-1;
            scheme(u,q, f, i, j, n, i-1, i-1, j+1, j-1, x, y, dtdx2, dtdy2, dt2, dt, b)
          
        #Update corners
        i = 0; j = 0
        scheme(u, q, f, i, j, n, i+1, i+1, j+1, j+1, x, y, dtdx2, dtdy2, dt2, dt, b)

        i = 0; j = Ny-1
        scheme(u, q, f, i, j, n, i+1, i+1, j-1, j-1, x, y, dtdx2, dtdy2, dt2, dt, b)

        i = Nx-1; j = 0
        scheme(u ,q, f, i, j, n, i-1, i-1, j+1, j+1, x, y, dtdx2, dtdy2, dt2, dt, b)

        i = Nx-1; j = Ny-1
        scheme(u, q, f, i, j, n, i-1, i-1, j-1, j-1, x, y, dtdx2, dtdy2, dt2, dt, b)
        
        
def scheme(u, q, f, i, j, n, i2, i3, j2, j3, x ,y, dtdx2, dtdy2, dt2, dt, b):
    """
    The numerical scheme. The same function can be used for all boundary conditions. By changing i2, i3..
    all different stencils can be recovered.
    
    Standard input is:
    i2 = i+1
    i3 = i-1

    j2 = j+1
    j3 = j-1

    These must be changed accordingly when working on the boundaries and corners
    """

    u[i,j,n+1] = 2*u[i,j,n] - (1 - 0.5*b*dt)*u[i,j,n-1] + \
    dtdx2*((q(x[i2],y[j]) + q(x[i],y[j]))*(u[i2,j,n] - u[i,j,n]) - (q(x[i],y[j]) + q(x[i3],y[j]))*(u[i,j,n] -u[i3,j,n])) + \
    dtdy2*((q(x[i],y[j2]) + q(x[i],y[j]))*(u[i,j2,n] - u[i,j,n]) - (q(x[i],y[j]) + q(x[i],y[j3]))*(u[i,j,n] -u[i,j3,n])) + \
    dt2*f(x[i],y[j],dt*n)
    
    u[i,j,n+1] /= 1 + 0.5*b*dt
    
def init_vec(u, I, V, x, y, dt):
    """
    Set initial conditions
    """
    u[:,:,0] = I(x,y)
    u[:,:,-1] = u[:,:,0] - dt*V(x,y)
            

    
def advance_vec(u, q, f, x, y, b, dx, dy, Nt, dt):
    """
    Advance the solution with time
    """
    dt2 = dt**2
    dtdx2 = dt2/(2*dx**2)
    dtdy2 = dt2/(2*dy**2)
    
    for n in xrange(0,Nt-1):
        scheme_vec(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt)

        #boundary conditions
        #j = 0;
        scheme_vec2(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   jstart = 0, jstop = 1, j2start = 1, j2stop = 2, j3start = 1, j3stop = 2)
        
        #j = Ny-1
        scheme_vec2(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   jstart = -1, jstop = None, j2start = -2, j2stop = -1, j3start = -2, j3stop = -1)

        
        #i = 0
        scheme_vec(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   istart = 0, istop = 1, i2start = 1, i2stop = 2, i3start = 1, i3stop = 2)
               
        #i = Nx-1
        scheme_vec(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   istart = -1, istop = None, i2start = -2, i2stop = -1, i3start = -2, i3stop = -1)
        
        #Create corners
        #i = 0; j = 0
        scheme_vec(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   istart = 0, istop = 1, i2start = 1, i2stop = 2, i3start = 1, i3stop = 2,\
                   jstart = 0, jstop = 1, j2start = 1, j2stop = 2, j3start = 1, j3stop = 2)
               
        #i = 0; j = Ny-1
        scheme_vec(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   istart = 0, istop = 1, i2start = 1, i2stop = 2, i3start = 1, i3stop = 2,\
                   jstart = -1, jstop = None, j2start = -2, j2stop = -1, j3start = -2, j3stop = -1)
                       
        #i = Nx-1; j = 0
        scheme_vec(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   istart = -1, istop = None, i2start = -2, i2stop = -1, i3start = -2, i3stop = -1,\
                   jstart = 0, jstop = 1, j2start = 1, j2stop = 2, j3start = 1, j3stop = 2)
                       
               
        #i = Nx-1; j = Ny-1
        scheme_vec(u, n, q, f, x, y, dtdx2, dtdy2, dt2, b, dt, \
                   istart = -1, istop = None, i2start = -2, i2stop = -1, i3start = -2, i3stop = -1,\
                   jstart = -1, jstop = None, j2start = -2, j2stop = -1, j3start = -2, j3stop = -1)
        
        
               

                   
                   
def scheme_vec(u, n, q, f, x ,y, dtdx2, dtdy2, dt2, b, dt, \
               istart = 1, istop = -1, i2start = 2, i2stop = None, i3start = 0, i3stop = -2,\
               jstart = 1, jstop = -1, j2start = 2, j2stop = None, j3start = 0, j3stop = -2):

    """
    The numerical scheme. We have impemented this version that  only makes it necessary to use two different functions for all boundary conditions. 

    The stencil can be changed by changing istart, istop... which are the values the different arrays go from-to. Thus this function works for almost all boundary conditions, except the horizontal ones, j =0 and j = Nx
    """
               
    u[istart:istop,jstart:jstop,n+1] = \
    2*u[istart:istop,jstart:jstop,n] - (1 - 0.5*b*dt)*u[istart:istop,jstart:jstop,n-1] + \
    dtdx2*((q(x[i2start:i2stop],y[jstart:jstop]) + q(x[istart:istop],y[jstart:jstop]))\
           *(u[i2start:i2stop,jstart:jstop,n] - u[istart:istop,jstart:jstop,n]) \
           - (q(x[istart:istop],y[jstart:jstop]) + q(x[i3start:i3stop],y[jstart:jstop]))\
           *(u[istart:istop,jstart:jstop,n] - u[i3start:i3stop,jstart:jstop,n])) + \
    dtdy2*((q(x[istart:istop],y[j2start:j2stop]) + q(x[istart:istop],y[jstart:jstop]))\
           *(u[istart:istop,j2start:j2stop,n] - u[istart:istop,jstart:jstop,n]) \
           - (q(x[istart:istop],y[jstart:jstop]) + q(x[istart:istop],y[j3start:j3stop]))\
           *(u[istart:istop,jstart:jstop,n] -u[istart:istop,j3start:j3stop,n])) + \
    dt2*f(x[istart:istop],y[jstart:jstop], dt*n)
    u[istart:istop,jstart:jstop,n+1] /=  1 + 0.5*b*dt
    

def scheme_vec2(u, n, q, f, x ,y, dtdx2, dtdy2, dt2, b, dt, \
               istart = 1, istop = -1, i2start = 2, i2stop = None, i3start = 0, i3stop = -2,\
               jstart = 1, jstop = -1, j2start = 2, j2stop = None, j3start = 0, j3stop = -2):

    """
    The numerical scheme. this function is made to make sure that the horizontal boundary conditions work, by reshaping the functions to fit into the horisontal boundary, j=0 and j = Nx in u[...]. 
    """
         
               
    u[istart:istop,jstart:jstop,n+1] = \
    2*u[istart:istop,jstart:jstop,n] - (1 - 0.5*b*dt)*u[istart:istop,jstart:jstop,n-1] + \
    dtdx2*((q(x[i2start:i2stop],y[jstart:jstop]) + q(x[istart:istop],y[jstart:jstop])).reshape(-1,1)\
           *(u[i2start:i2stop,jstart:jstop,n] - u[istart:istop,jstart:jstop,n]) \
           - (q(x[istart:istop],y[jstart:jstop]) + q(x[i3start:i3stop],y[jstart:jstop])).reshape(-1,1)\
           *(u[istart:istop,jstart:jstop,n] - u[i3start:i3stop,jstart:jstop,n])) + \
    dtdy2*((q(x[istart:istop],y[j2start:j2stop]) + q(x[istart:istop],y[jstart:jstop])).reshape(-1,1)\
           *(u[istart:istop,j2start:j2stop,n] - u[istart:istop,jstart:jstop,n]) \
            - (q(x[istart:istop],y[jstart:jstop]) + q(x[istart:istop],y[j3start:j3stop])).reshape(-1,1)\
            *(u[istart:istop,jstart:jstop,n] -u[istart:istop,j3start:j3stop,n])) + \
            dt2*(f(x[istart:istop],y[jstart:jstop], dt*n)).reshape(-1,1)
    u[istart:istop,jstart:jstop,n+1] /=  1 + 0.5*b*dt

    
    
if __name__ == '__main__':
    print "something"

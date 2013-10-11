import nose.tools as nt
import pylab as p


def solver(b, Lx, dx, Ly, dy, T, dt, version="scalar"):
    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    u = p.zeros((Nx,Ny,Nt))

    x = p.linspace(0, Lx, Nx)  # mesh points in x dir
    y = p.linspace(0, Ly, Ny)  # mesh points in y dir

    #Create initial conditions stuff
    if version == 'vectorized':
        init_vec(u, x, y, dt)
        advance_vec(u, x, y, b, dx, dy, Nt, dt)
    else:
        init(u, Nx, Ny, x, y, dt)
        advance(u, x, y, b, Nx, dx, Ny, dy, Nt ,dt)
                
    return u


    
    
def init(u, Nx, Ny, x, y, dt):

    for i in xrange(1, Nx-1):
        for j in xrange(1, Ny-1):
            u[i,j,0] = I(x[i],y[j])

    for i in xrange(1, Nx-1):
        for j in xrange(1, Ny-1):
            u[i,j,-1] = u[i,j,0] - dt*V(x[i],y[j])
            


def advance(u, x, y, b, Nx, dx, Ny, dy, Nt, dt):
    dt2 = dt**2
    dtdx2 = dt2/(2*dx**2)
    dtdy2 = dt2/(2*dy**2)
    
    for n in xrange(0,Nt-1):
        for i in xrange(1, Nx-1):
            for j in xrange(1, Ny-1):
                scheme(u,i,j,n, i+1, i-1, j+1, j-1, x, y, dtdx2, dtdy2, dt2)

        #Update boundaries
        for i in xrange(1,Nx-1):
            j = 0;
            scheme(u,i,j,n, i+1, i-1, j+1, j+1, x, y, dtdx2, dtdy2, dt2)
            
            j = Ny-1;
            scheme(u,i,j,n, i+1, i-1, j-1, j-1, x, y, dtdx2, dtdy2, dt2)
            
        for j in xrange(1,Ny-1):
            i = 0;
            scheme(u,i,j,n, i+1, i+1, j+1, j-1, x, y, dtdx2, dtdy2, dt2)
            
            i = Nx-1;
            scheme(u,i,j,n, i-1, i-1, j+1, j-1, x, y, dtdx2, dtdy2, dt2)
            
            #return u 

        #Update corners
        i = 0; j = 0
        scheme(u,i,j,n, i+1, i+1, j+1, j+1, x, y, dtdx2, dtdy2, dt2)

        i = 0; j = Ny-1
        scheme(u,i,j,n, i+1, i+1, j-1, j-1, x, y, dtdx2, dtdy2, dt2)

        i = Nx-1; j = 0
        scheme(u,i,j,n, i-1, i-1, j+1, j+1, x, y, dtdx2, dtdy2, dt2)

        i = Nx-1; j = Ny-1
        scheme(u,i,j,n, i-1, i-1, j-1, j-1, x, y, dtdx2, dtdy2, dt2)

        
def scheme(u, i, j, n, i2, i3, j2, j3, x ,y, dtdx2, dtdy2, dt2):
    """
    Standard input is:
    i2 = i+1
    i3 = i-1

    j2 = j+1
    j3 = j-1
    """
    
    """
    u[i,j,n+1] = 2*u[i,j,n] - (1 - b*dt)*u[i,j,n-1] + \\
    dydx2*((q(i+1,j) + q(i,j))*(u[i+1,j,n] - u(i,j,n)) - (q(i,j) + q(i-1,j)*(u[i,j,n] -u[i-1,j,n]))) + \\
    dtdy2*((q(i,j+1) + q(i,j))*(u[i,j+1,n] - u[i,j,n]) - (q(i,j) + q(i,j-1)*(u[i,j,n] -u[i,j-1,n]))) + \\
    dt2*f(i,j,n)
    u[i,j,n+1] /= 1 + b*dt
    """
    u[i,j,n+1] = 2*u[i,j,n] - (1 - 0.5*b*dt)*u[i,j,n-1] + \
    dtdx2*((q(x[i2],y[j]) + q(x[i],y[j]))*(u[i2,j,n] - u[i,j,n]) - (q(x[i],y[j]) + q(x[i3],y[j])*(u[i,j,n] -u[i3,j,n]))) + \
    dtdy2*((q(x[i],y[j2]) + q(x[i],y[j]))*(u[i,j2,n] - u[i,j,n]) - (q(x[i],y[j]) + q(x[i],y[j3])*(u[i,j,n] -u[i,j3,n]))) + \
    dt2*f(x[i],y[j],n)
    u[i,j,n+1] /= 1 + 0.5*b*dt
    

def init_vec(u, x, y, dt):
    u[:,:,0] = I(x,y)
    u[:,:,-1] = u[:,:,0] - dt*V(x,y)
            

    
def advance_vec(u, x, y, b, dx, dy, Nt, dt):
    dt2 = dt**2
    dtdx2 = dt2/(2*dx**2)
    dtdy2 = dt2/(2*dy**2)
    
    for n in xrange(0,Nt-1):
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b)

        #boundary conditions
        #j = 0;
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b,\
                   jstart = 0, jstop = 1, j2start = 1, j2stop = 2, j3start = 1, j3stop = 2)
    
        #j = Ny-1
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b,\
                   jstart = -1, jstop = None, j2start = -2, j2stop = -1, j3start = -2, j3stop = -1)

        #i = 0
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b,\
                   istart = 0, istop = 1, i2start = 1, i2stop = 2, i3start = 1, i3stop = 2)
                   
        #i = Nx-1
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b,\
                   istart = -1, istop = None, i2start = -2, i2stop = -1, i3start = -2, i3stop = -1)
        
        #Create corners
        #i = 0; j = 0
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b, \
                   istart = 0, istop = 1, i2start = 1, i2stop = 2, i3start = 1, i3stop = 2,\
                   jstart = 0, jstop = 1, j2start = 1, j2stop = 2, j3start = 1, j3stop = 2)
               
        #i = 0; j = Ny-1
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b, \
                   istart = 0, istop = 1, i2start = 1, i2stop = 2, i3start = 1, i3stop = 2,\
                   jstart = -1, jstop = None, j2start = -2, j2stop = -1, j3start = -2, j3stop = -1)
                       
        #i = Nx-1; j = 0
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b, \
                   istart = -1, istop = None, i2start = -2, i2stop = -1, i3start = -2, i3stop = -1,\
                   jstart = 0, jstop = 1, j2start = 1, j2stop = 2, j3start = 1, j3stop = 2)
                       
               
        #i = Nx-1; j = Ny-1
        scheme_vec(u, n, x, y, dtdx2, dtdy2, dt2, b, \
                   istart = -1, istop = None, i2start = -2, i2stop = -1, i3start = -2, i3stop = -1,\
                   jstart = -1, jstop = None, j2start = -2, j2stop = -1, j3start = -2, j3stop = -1)

        
               

                   
                   
def scheme_vec(u, n, x ,y, dtdx2, dtdy2, dt2, b, \
               istart = 1, istop = -1, i2start = 2, i2stop = None, i3start = 0, i3stop = -2,\
               jstart = 1, jstop = -1, j2start = 2, j2stop = None, j3start = 0, j3stop = -2):

    """
        print u[istart:istop,jstart:jstop,n+1].shape

    print  (2*u[istart:istop,jstart:jstop,n] - (1 - 0.5*b*dt)*u[istart:istop,jstart:jstop,n-1] + \
    dtdx2*((q(x[i2start:i2stop],y[jstart:jstop]) + q(x[istart:istop],y[jstart:jstop]))\
           *(u[i2start:i2stop,jstart:jstop,n] - u[istart:istop,jstart:jstop,n]) \
           - (q(x[istart:istop],y[jstart:jstop]) + q(x[i3start:i3stop],y[jstart:jstop])\
              *(u[istart:istop,jstart:jstop,n] - u[i3start:i3stop,jstart:jstop,n]))) + \
    dtdy2*((q(x[istart:istop],y[j2start:j2stop]) + q(x[istart:istop],y[jstart:jstop]))\
           *(u[istart:istop,j2start:j2stop,n] - u[istart:istop,jstart:jstop,n]) \
    - (q(x[istart:istop],y[jstart:jstop]) + q(x[istart:istop],y[j3start:j3stop])\
           *(u[istart:istop,jstart:jstop,n] -u[istart:istop,j3start:j3stop,n]))) + \
        dt2*f(x[istart:istop],y[jstart:jstop], dt*n)).shape
        """

    u[istart:istop,jstart:jstop,n+1] = \
          2*u[istart:istop,jstart:jstop,n] - (1 - 0.5*b*dt)*u[istart:istop,jstart:jstop,n-1] + \
    dtdx2*((q(x[i2start:i2stop],y[jstart:jstop]) + q(x[istart:istop],y[jstart:jstop]))\
           *(u[i2start:i2stop,jstart:jstop,n] - u[istart:istop,jstart:jstop,n]) \
           - (q(x[istart:istop],y[jstart:jstop]) + q(x[i3start:i3stop],y[jstart:jstop])\
              *(u[istart:istop,jstart:jstop,n] - u[i3start:i3stop,jstart:jstop,n]))) + \
    dtdy2*((q(x[istart:istop],y[j2start:j2stop]) + q(x[istart:istop],y[jstart:jstop]))\
           *(u[istart:istop,j2start:j2stop,n] - u[istart:istop,jstart:jstop,n]) \
    - (q(x[istart:istop],y[jstart:jstop]) + q(x[istart:istop],y[j3start:j3stop])\
           *(u[istart:istop,jstart:jstop,n] -u[istart:istop,j3start:j3stop,n]))) + \
        dt2*f(x[istart:istop],y[jstart:jstop], dt*n)
        
    u[istart:istop,jstart:jstop,n+1] /= 1 + 0.5*b*dt
    
            
    """
    u[i:j,i:j,n+1] = 2*u[i:j,1:-2,n] - (1 - 0.5*b*dt)*u[i:j,i:j,,n-1] + \
    dydx2*((q(x[i:j],y[i:j]) + q(x[i:j],y[i:j]))*(u[i2:j2,i:j,n] - u[i:j,i:j,n]) - (q(x[i:j],y[i:j]) + q(x[i3:j3],y[i3:j3])*(u[i:j,i:j,n] -u[i3:j3,i:j,n]))) + \
    dtdy2*((q(x[i:j],y[i2:j2]) + q(x[i:j],y[i:j]))*(u[i:j,i2:j2,n] - u[i:j,i:j,n]) - (q(x[i:j],y[i:j]) + q(x[i:j],y[i3:j3])*(u[i:j,i:j,n] -u[i:j,i3:j3,n]))) + \
    dt2*f(x[i:j],y[i:j],n)
    u[i:j,i:j,n+1] /= 1 + 0.5*b*dt
"""    

    
    
def V(x,y):
    return 7

def I(x, y):
    return 3

def q(x, y):
    return 15000

def f(x, y, n):
    return 0



def testConstant():
    b = 0.1
    Lx = 10.
    Ly = 10.

    dx = 0.1
    dy = 0.1
    dt = 0.1
    T = 5.
    
    u = solver(b, Lx, dx, Ly, dy, T, dt, version="vectorize")
    


if __name__ == '__main__':
    print "something"

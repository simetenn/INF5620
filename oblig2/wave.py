import nose.tools as nt
import pylab as p


def solver(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="scalar"):
    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))


    if version == 'vectorized':
        advance = advance_vec
    else version == 'scalar':
        advance = advance

    
    u = p.zeros(Nx,Ny,Nt)

    x = linspace(0, Lx, Nx+1)  # mesh points in x dir
    y = linspace(0, Ly, Ny+1)  # mesh points in y dir

    #Create initial conditions stuff
    
    
    
def init(u, I, V, Nx, Ny, x, y, q, f, b, dx, dy, dt, scheme):

    for i in xrange(1, Nx-1):
        for j in xrange(1, Ny-1):
            u[i,j,0] = I(x[i],y[j])

    for i in xrange(1, Nx-1):
        for j in xrange(1, Ny-1):
            
            


def advance(u, Nx, Ny, x, y, q, f, b, dx, dy, dt, scheme):
    dt2 = dt**2
    dtdx2 = dt2/(2*dx**2)
    dtdy2 = dt2/(2*dy**2)
    
    for n in xrange(1,Nt-1):
        for i in xrange(1, Nx-1):
            for j in xrange(1, Ny-1):
                scheme(u,i,j,n, i+1, i-1, j+1, j-1, x, y)

        #Oppdatere kanter
        for i in xrange(1,Nx-1):
            j = 0;
            scheme(u,i,j,n, i+1, i-1, j+1, j+1, x, y)
            
            j = Nx-1;
            scheme(u,i,j,n, i+1, i-1, j-1, j-1, x, y)
            
        for j in xrange(1,Ny-1):
            i = 0;
            scheme(u,i,j,n, i+1, i+1, j+1, j-1, x, y)
            
            i = Nx-1;
            scheme(u,i,j,n, i-1, i-1, j+1, j-1, x, y)
            
            #return u 


def scheme(u, i, j, n, i2, i3, j2, j3,x ,y):
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
    u[i,j,n+1] = 2*u[i,j,n] - (1 - b*dt)*u[i,j,n-1] + \\
    dydx2*((q(x[i2],y[j]) + q(x[i],y[j]))*(u[i2,j,n] - u[i,j,n]) - (q(x[i],y[j]) + q(x[i3],y[j])*(u[i,j,n] -u[i3,j,n]))) + \\
    dtdy2*((q(x[i],y[j2]) + q(x[i],y[j]))*(u[i,j2,n] - u[i,j,n]) - (q(x[i],y[j]) + q(x[i],y[j3])*(u[i,j,n] -u[i,j3,n]))) + \\
    dt2*f(x[i],y[j],n)
    u[i,j,n+1] /= 1 + b*dt
    

#def advance_vec():

#
#def init_vec(u, I, V, Nx, Ny, x, y, q, f, b, dx, dy, dt, scheme):
#    u[:,:,0] = I(x,y)
    

    
def V(x,y):
    return 0

def I(x,y):
    return 0

def q(x,y):
    return 0



if __name__ == '__main__':
    print "something"

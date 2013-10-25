from wave import *
from pulse import pulse


def test_constant_solution():
    b = 0
    Lx = 1.
    Ly = 1.

    dx = 0.1
    dy = 0.1
    dt = 0.1
    T = 1.

    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    x = p.linspace(0, Lx, Nx)
    y = p.linspace(0, Ly, Ny)
    t = p.linspace(0, T, Nt)
    
    def V(x, y):
        return 0

    def I(x, y):
        return 10

    def q(x, y):
        return 5

    def f(x, y, n):
        return 0

    
    u = solver(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="scalar")

    def exact_solution(x, y, t):
         return 10
    
    u_e = exact_solution(x, y, t)
    difference = abs(u_e - u).max()
    print "Largest difference: ", difference
    nt.assert_almost_equal(difference, 0, places=15)



def test_constant_solution_vec():
    b = 0
    Lx = 10.
    Ly = 10.

    dx = 0.1
    dy = 0.1
    dt = 0.01
    T = 1.

    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    x = p.linspace(0, Lx, Nx)
    y = p.linspace(0, Ly, Ny)
    t = p.linspace(0, T, Nt)
    
    def V(x, y):
        return 0

    def I(x, y):
        return 10

    def q(x, y):
        return 5

    def f(x, y, n):
        return 0

    
    u = solver(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")

    
    
    def exact_solution(x, y, t):
         return 10

    u_e = exact_solution(x, y, t)
    difference = abs(u_e - u).max()
    print "Largest difference: ", difference
    nt.assert_almost_equal(difference, 0, places=15)


def test_constant_solution_vec():
    b = 0
    Lx = 10.
    Ly = 10.

    dx = 0.1
    dy = 0.1
    dt = 0.01
    T = 1.

    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    x = p.linspace(0, Lx, Nx)
    y = p.linspace(0, Ly, Ny)
    t = p.linspace(0, T, Nt)
    
    def V(x, y):
        return 0

    def I(x, y):
        return 10

    def q(x, y):
        return 5

    def f(x, y, n):
        return 0



    u = solver(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")

    
    def exact_solution(x, y, t):
         return 10
    
    u_e = exact_solution(x, y, t)
    difference = abs(u_e - u).max()
    print "Largest difference: ", difference
    nt.assert_almost_equal(difference, 0, places=15)



def test_plug():
    """Check that an initial plug is correct back after one period."""
  
    b = 0
    Lx = 1.
    Ly = 1.

    dx = 0.05
    dy = 0.05
    dt = 0.1
    T = 4.

    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    x = p.linspace(0, Lx, Nx)
    y = p.linspace(0, Ly, Ny)
    t = p.linspace(0, T, Nt)
    
    def V(x, y):
        return 0

    def Ix(x, y):
        result = p.zeros((len(x),len(y)))
        result[:,:] = 1
        result[p.where(abs(x-Lx/2.0) > 0.1),:] = 0
        return result

    def Isx(x,y):
        if abs(x-Lx/2.0) > 0.1:
            return 0
        else:
            return 1


    def Iy(x, y):
        result = p.zeros((len(x),len(y)))
        result[:,:] = 1
        result[:,p.where(abs(y-Ly/2.0) > 0.1)] = 0
        return result

    def Isy(x,y):
        if abs(y-Ly/2.0) > 0.1:
            return 0
        else:
            return 1
        

    def q(x, y):
        return 5

    def f(x, y, n):
        return 0


    u_vx = solver(Ix, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")
    u_sx = solver(Isx, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="scalar")

    u_vy = solver(Iy, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")
    u_sy = solver(Isy, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="scalar")
    
 
    
    diff = abs(u_sx - u_vx).max()
    nt.assert_almost_equal(diff, 0, places=13)
    diff = abs(u_sy - u_vy).max()
    nt.assert_almost_equal(diff, 0, places=13)






def test_undampened():
    b = 0
    Lx = 10.
    Ly = 10.
    T = 10.
 
    def V(x, y):
        return 0

    def I(x, y):
        return u_e(x, y, 0)
      
    def q(x, y):
        return 10

    def f(x, y, n):
        return 0


    def u_e(x,y,t):
        A = 8
        mx = .1
        my = .1
        kx = mx*p.pi/Lx
        ky = my*p.pi/Ly
        omega = 0.1
        x,y = p.meshgrid(x,y)
        return A*p.cos(kx*x)*p.cos(ky*y)*p.cos(omega*t)
            
    def u_e2(x,y,t):
        A = 1
        mx = 0.1
        my = 0.1
        kx = mx*p.pi/Lx
        ky = my*p.pi/Ly
        omega = 1

        return A*p.cos(kx*x)*p.cos(ky*y)*p.cos(omega*t)
 
        
        
    expected_rate = 2
    n = 10
    E = []
    c = 0.1
    h_list = p.linspace(5, 0.5, n)
    #h_list = [0.01]
    for h in h_list: 

        Nx = int(round(Lx/float(h)))
        Ny = int(round(Ly/float(h)))
        #Nt = int(round(T/float(c*h)))
        
        
        
        x = p.linspace(0, Lx, Nx)
        y = p.linspace(0, Ly, Ny)
        #t = p.linspace(0, T, Nt)
        
        u = solver(I, V, q, f, b, Lx, h, Ly, h, T, c*h, version="vectorized")
        #xx,yy = p.meshgrid(x,y)
        v_e = u_e(x,y,T)

        
        # print
        # print v_e[:,:,-1]
        #print "----------------"
        #print u[:,:,-1]
        #sys.exit()
        
        E.append(abs(v_e - u[:,:,-1]).max())

        #xx,yy = p.meshgrid(x,y)
        #p.np.save("u", u)
    #p.np.save("h",B1(xx,yy))
    #p.np.save("x",x)
    #p.np.save("y",y)
    E = p.array(E)
    rate = p.zeros(n-1)
    for i in xrange(1, n):
        rate[i-1] = p.log(E[i-1]/E[i])/p.log(h_list[i-1]/h_list[i])

    print E/h_list**2
    diff = abs(expected_rate - rate[-1])
    nt.assert_almost_equal(diff, 0 ,places=1)
 






def test_dampened():
    b = 1
    Lx = 10.
    Ly = 10.
    T = 10.
 
    
    def V(x, y):
        return 0

    def I(x, y):
        return u_e(x, y, 0)
    
    def q(x, y):
        return 10

    def f(x, y, n):
        return 0

    def u_e(x,y,t):
        A = 1
        B = 0
        mx = 1
        my = 1
        kx = mx*p.pi/Lx
        ky = my*p.pi/Ly
        c = b/2. 
        omega = p.sqrt(kx**2*q(x,y) + ky**2*q(x,y) - c**2)
        x,y = p.meshgrid(x,y)
        return (A*p.cos(omega*t) + B*p.sin(omega*t))*p.cos(kx*x)*p.cos(ky*y)*p.exp(-c*t)

        
    expected_rate = 2
    n = 20
    E = []
    c = 0.1
    h_list = p.linspace(5, 0.05, n)
    #h_list = [0.01]
    for h in h_list: 
        Nx = int(round(Lx/float(h)))
        Ny = int(round(Ly/float(h)))

        x = p.linspace(0, Lx, Nx)
        y = p.linspace(0, Ly, Ny)
        
        u = solver(I, V, q, f, b, Lx, h, Ly, h, T, c*h, version="vectorized")
        v_e = u_e(x,y,T)
        
        E.append(abs(v_e - u[:,:,-1]).max())

    E = p.array(E)
    rate = p.zeros(n-1)
    for i in xrange(1, n):
        rate[i-1] = p.log(E[i-1]/E[i])/p.log(h_list[i-1]/h_list[i])

        #print E/h_list**2
    print rate
    diff = abs(expected_rate - rate[-1])
    nt.assert_almost_equal(diff, 0 ,places=1)
 


def test_mms():
    b = 1
    Lx = 10.
    Ly = 10.
    T = 10.
 
    
    def V(x, y):
        return 0

    def I(x, y):
        return u_e(x, y, 0)
    
    def q(x, y):
        return 10

    def f(x, y, n):
        return 0

    def u_e(x,y,t):
        A = 1
        B = 0
        mx = 1
        my = 1
        kx = mx*p.pi/Lx
        ky = my*p.pi/Ly
        c = b/2. 
        omega = p.sqrt(kx**2*q(x,y) + ky**2*q(x,y) - c**2)
        x,y = p.meshgrid(x,y)
        return (A*p.cos(omega*t) + B*p.sin(omega*t))*p.cos(kx*x)*p.cos(ky*y)*p.exp(-c*t)

        
    expected_rate = 2
    n = 20
    E = []
    c = 0.1
    h_list = p.linspace(5, 0.05, n)
    #h_list = [0.01]
    for h in h_list: 
        Nx = int(round(Lx/float(h)))
        Ny = int(round(Ly/float(h)))

        x = p.linspace(0, Lx, Nx)
        y = p.linspace(0, Ly, Ny)
        
        u = solver(I, V, q, f, b, Lx, h, Ly, h, T, c*h, version="vectorized")
        v_e = u_e(x,y,T)
        
        E.append(abs(v_e - u[:,:,-1]).max())

    E = p.array(E)
    rate = p.zeros(n-1)
    for i in xrange(1, n):
        rate[i-1] = p.log(E[i-1]/E[i])/p.log(h_list[i-1]/h_list[i])

    print rate
    diff = abs(expected_rate - rate[-1])
    nt.assert_almost_equal(diff, 0 ,places=1)




    

def physical():
    b = 1
    Lx = 2
    Ly = 2

    dx = 0.1
    dy = 0.1
    dt = 0.01
    T = 1

    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    x = p.linspace(0, Lx, Nx)
    y = p.linspace(0, Ly, Ny)
    t = p.linspace(0, T, Nt)
    
    def V(x, y):
        return 0

    
    def I(x, y):
        I0 = 2
        Ia = 1
        Im = 1
        Is = 0.5*p.sqrt(2)
        return I0 + Ia*p.exp(-((x - Im)/Is)**2) 


    def I2(x, y):
        I0 = 2
        Ia = 1
        Im = 1
        Is = 0.5
        return I0 + Ia*p.exp(-((x - Im)/Is)**2-((y.reshape(-1,1) - Im)/Is)**2) 


    def B1(x,y):
        B0 = 0
        Ba = 1
        Bmy = 1
        Bmx = 1
        Bs = p.sqrt(2)*0.5
        b = 1
        return B0 + Ba*p.exp(-((x - Bmx)/Bs)**2-((y - Bmy)/(b*Bs))**2) 


    def B2(x,y):
        B0 = 0
        Ba = 1
        Bmy = 5
        Bmx = 5
        Bs = 3
        b = 1
        
        index = 0 <  p.sqrt(x**2 + y**2)
        index2 = p.sqrt(x**2 + y**2) <= Bs
        index = index+index2
        
        results =  B0 + Ba*p.cos(p.pi*((x - Bmx)/(2*Bs)))*p.cos(p.pi*(y - Bmy)/(2*Bs)) 

        
        for i in xrange(len(x)):
            for j in xrange(len(y)):
                if (0 > p.sqrt((x[i])**2 + (y[j])**2) >= Bs):
                    results[i,j] = B0
        
        """
        results = p.zeros((len(x),len(y)))
        index = x > Bmx - Bs 
        index2 = x < (Bmx + Bs)
        xindex = p.invert(index - index2)
        index = y > Bmy - Bs 
        index2 = y < (Bmy + Bs)
        yindex = p.invert(index - index2)
        results[:,:] =  B0
        results[xindex*yindex] = B0 + Ba
        """
        #print index
        #print index2
        #results[not (index[0])] = B0
        #results[not index2] = B0
        #print results
        return results
    
    def B3(x,y):
        B0 = 0
        Ba = 1
        Bmy = 5
        Bmx = 5
        Bs = 1
        b = 1

        results = p.zeros((len(x),len(y)))
        index = x > Bmx - Bs 
        index2 = x < (Bmx + Bs)
        xindex = p.invert(index - index2)
        index = y > Bmy -b*Bs 
        index2 = y < (Bmy + b*Bs)
        yindex = p.invert(index - index2)
        results[:,:] =  B0
        results[xindex*yindex] = B0 + Ba
        return results
    


    
    def H(x,y):
        I0 = 2
        return I0 - B1(x,y)

    def q(x, y):
        return 9.81*H(x,y)


    def f(x, y, n):
        return 0



    u = solver(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")

    #print u[:,:,-1]
    #u = solver(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vecto")

    #print u[:,:,-1]
    #print y.reshape(-1,1)
    xx,yy = p.meshgrid(x,y)
    p.np.save("u", u)
    p.np.save("h",B1(x,y.reshape(-1,1)))
    p.np.save("x",x)
    p.np.save("y",y)
    
    
    


    
if __name__ == '__main__':
    #test_constant_solution_vec()
    physical()
    #test_plug()
    #test_constant_solution()
    #test_dampened()

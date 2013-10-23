from wave import *
    
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

if __name__ == '__main__':
    test_constant_solution_vec()
    test_constant_solution()

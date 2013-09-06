import nose.tools as nt
import pylab as p
from math import log

class Solver(object):
    """
    A solver class for solving differential equations
    Contains the numerical data, such as z, v ...
    """

    def __init__(self, problem):
        """
        Initializing the problem, using the Problem class
        """
        self.problem = problem


    def forces(self):
        """
        Calculating the forces, must be done after solve() have been run
        """
        def Fg():
            """
            Gravity force
            """
            return -self.problem.m*self.problem.g*p.ones(len(self.t))
        
        def Fb():
            """
            Buoyancy force
            """
            return self.problem.rho*self.problem.g*self.problem.V*p.ones(len(self.t))

        def Fd():
            """
            Drag force
            """
            return -0.5*self.problem.Cd(self.t)*self.problem.rho*self.problem.A(self.t)*abs(self.v)*self.v

        self.Fg_list  = Fg()
        self.Fd_list  = Fd()
        self.Fb_list  = Fb()
        
        
    def solve(self):
        """
        Solve the physical problem given
        Using a Crank-Nicolson scheme for the velocity and a Euler-Chromer scheme to update the position
        """
        
        self.z = p.zeros(self.problem.N + 1)
        self.z[0] = self.problem.z0 
        self.v = p.zeros(self.problem.N + 1)
        self.v[0] = self.problem.v0 
        self.t = p.linspace(0, self.problem.T, self.problem.N + 1) 
        
        dt = self.problem.dt

        def G(f, n):
            """
            Generalized arithmetic mean of a function f
            """
            return 0.5*(f(self.t[n + 1]) + f(self.t[n]))



        for n in xrange(self.problem.N):
            self.v[n+1] = (self.v[n] + G(self.problem.b, n)*dt + G(self.problem.c, n)*dt*self.problem.F(n))/(1 + G(self.problem.a, n)*dt*abs(self.v[n]))
            self.z[n+1] = self.v[n+1]*dt + self.z[n]  


    def __call__(self):
        return self.v, self.t
        

        
        
    def plot(self, elements):
        """
        Plot the results from the solver
        Input given as a list, the content gives which functions to plot

        0 => Position
        1 => Velocity
        2 => Drag force
        3 => Gravity force
        4 => Buoyancy force

        Saves the plot in a file with the appropriate name
        """
        
        for element in elements:
            if element == 0:
                p.figure()
                p.plot(self.t, self.z)
                p.ylabel("Position, m")
                p.xlabel("time, t")
                p.savefig("position.png")
                
            if element == 1:
                p.figure()
                p.plot(self.t, self.v)
                p.ylabel("Velocity, m/s")
                p.xlabel("time, t")
                p.savefig("velocity.png")

            if element == 2:
                p.figure()
                p.plot(self.t, self.Fd_list)
                p.ylabel("Drag, N")
                p.xlabel("time, t")
                p.savefig("drag.png")
                
            if element == 3:
                p.figure()
                p.plot(self.t, self.Fg_list)
                p.ylabel("Gravity, N")
                p.xlabel("time, t")  
                p.savefig("gravity.png")
                
            if element == 4:
                p.figure()
                p.plot(self.t, self.Fb_list)
                p.ylabel("Buoyancy, N")
                p.xlabel("time, t")
                p.savefig("buoyancy.png")

                
class Problem(object):
    """
    Class containing the physical data for a given problem
    """
    
    def __init__(self, dt, T, v0, z0, rho, m, V, test=False, tp = False):
        """
        Initialize the problem
        """
        
        self.initialize(dt, T, v0, z0, rho, m, V, test, tp)
    

    def initialize(self, dt, T, v0, z0, rho, m, V, test=False, tp = False):
        """
        Initialize the problem
        """
        self.dt = float(dt)
        self.T = T
        self.tp = tp
        self.v0 = v0
        self.N = int(round(T/self.dt)) 
        self.v0 = v0
        self.z0 = z0
        self.rho = rho
        self.m = m
        self.V = V
        self.g = 9.81                 #m/s^2

        self.test = test
        
        self.alpha = 1.5
        self.beta = 0


    def Cd(self, t):
        """
        The drag coefficient of a skydiver, changes when the parachute is deployed, at t=tp
        If tp = False, it returns the drag coefficient of a skydiver withouth a parachute
        Works for both arrays and single numbers
        """
        if not isinstance(t, p.ndarray):
            t = p.array([t],float)

        a = t.copy()
        if self.tp:
            a[t < self.tp] = 1.4
            a[t >= self.tp] = 1.8
        else:
            a[:] = 1.4
        return a

        
    def A(self, t):
        """
        The area of the skydiver, m^2. Changes when the parachute is deployed, at t=tp
        If tp = False, it returns the area of a skydiver withouth a parachute
        Works for both arrays and single numbers
        """
        if not isinstance(t, p.ndarray):
            t = p.array([t],float)

        a = t.copy()
        if self.tp:
            a[t < self.tp] = 0.5
            a[t >= self.tp] = 44.
        else:
            a[:] = 0.5
        return a

    
    def a(self, t):
        """
        The a function for the "constant" a, that varies with time  
        """
        return 0.5*(self.Cd(t)*self.rho*self.A(t))/self.m

    
    def b(self, t):
        """
        The function for the constant b
        """
        return self.g*(self.rho*self.V/self.m - 1)

    
    def c(self, t):
        """
        The function for the constant c
        """
        return 1./self.m

    
    def F(self, n):
        """
        A source term, containing an MMS

        "sin"  => v(t) = sin(alpha*t) 
        "ax+b" => v(t) = alpha*t + beta
        anything else => the actual solution, the source term is set to 0
        """
        if self.test == "ax+b": 
            return (self.alpha + self.a(n*self.dt)*self.alpha*n*(n+1)*self.dt**2*abs(self.alpha) - self.b(n*self.dt))/self.c(n*self.dt) 
        elif self.test == "sin":
            return (p.cos(self.alpha*n*self.dt) + self.a(n*self.dt)*abs(p.sin(self.alpha*n*self.dt))*p.sin(self.alpha*(n+1)*self.dt) - self.b(n*self.dt))/self.c(n*self.dt)
        return 0
    



def test_linear_solution():
    """
    Test if a MMS linear solution gives correct answe using the Solver()
    """
    dt = 0.05         #Time step, s
    T = 30.           #End time for simulation, s
    v0 = 0.           #Initial velocity, m/s
    z0 = 3000.        #Initial position, m
    rho = 1.          #Density of air, kg/m^3
    m = 100.          #Mass of skydiver, kg
    V = 0.07          #Volume of skydiver, kg

    
    #Solve the problem with a linear MMS
    P = Problem(dt, T, v0, z0, rho, m, V, test="ax+b")
    s = Solver(P)
    s.solve()
    v, t = s()

    
    def exact_solution(t):
        return P.alpha*t + P.beta

    
    v_e = exact_solution(t)
    difference = abs(v_e - v).max()
    print "Largest difference: ", difference
    nt.assert_almost_equal(difference, 0, places=12)


    

def test_convergence_rate():
    """
    Test if the convergence rate for the solver is reasonable
    """
    T = 30.           #End time for simulation, s
    v0 = 0.           #Initial velocity, m/s
    z0 = 3000.        #Initial position, m
    rho = 1.          #Density of air, kg/m^3
    m = 100.          #Mass of skydiver, kg
    V = 0.07          #Volume of skydiver, kg
    n = 11            #Nr dt steps
    expected_rate = 2 #Expected convergence rate for the problem
    tol = 0.1         #Tolerance

    #Initialize an instance of the class
    P = Problem(1, T, v0, z0, rho, m, V)

    def exact_solution(t):
        return p.sin(P.alpha*t)

    
    E = []
    dt_list = p.linspace(0.5, 0.005, n)
    for dt in dt_list: 
        #Solve the problem with a sine MMS
        P.initialize(dt, T, v0, z0, rho, m, V, test = "sin")
        s = Solver(P)
        s.solve()
        v, t = s()
        v_e = exact_solution(t)
        E.append(abs(v_e - v).max())
    
    rate = p.zeros(n-1)
    
    for i in xrange(1, n):
        rate[i-1] = p.log(E[i-1]/E[i])/p.log(dt_list[i-1]/dt_list[i])

    diff = abs(expected_rate - rate[-1])
    nt.assert_almost_equal(diff, 0 ,places=1)
        

    

if __name__ == "__main__":
    dt = 0.05         #Time step, s
    T = 300.          #End time for simulation, s
    v0 = 0.           #Initial velocity, m/s
    z0 = 3000.        #Initial position, m
    rho = 1.          #Density of air, kg/m^3
    m = 100.          #Mass of skydiver, kg
    V = 0.07          #Volume of skydiver, kg
    
    P = Problem(dt, T, v0, z0, rho, m, V, tp = 60.)
    S = Solver(P)
    S.solve()
    S.forces()
    S.plot([0,1,2,3,4])
        
    #test_convergence_rate()
    #test_linear_solution()


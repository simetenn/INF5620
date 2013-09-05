class Solver:
    """
    A class for solving differential equations
    """

    def __init__(self):
        """
        Initializing
        """


        
    def solve(self,function):
        """
        Solve the physical problem given
        """
        
        def problem(self,):
        """
        Initialize the problem
        """
        
    def plot(self, elements):
        """
        Plot the results from the solver
        """
        import pylab as p

        for element in elements:
            if element == 0:
                p.figure()
                p.plot(self.x,self.t)
                p.xlabel("Position, m/s")
                p.ylable("time, t")

            if element == 1:
                p.figure()
                p.plot(self.v,self.t)
                p.xlabel("Velocity, m/s")
                p.ylable("time, t")
        

            if element == 2:
                p.figure()
                p.plot(self.Fd,self.t)
                p.xlabel("Drag, N")
                p.ylable("time, t")

            if element == 3:
                    p.figure()
                    p.plot(self.Fg, self.t)
                    p.xlabel("Gravity, N")
                    p.ylable("time, t")
        
            if element == 4:
                p.figure()
                p.plot(self.Fg, self.t)
                p.xlabel("Buoyancy, N")
                p.ylable("time, t")


        
if name == "__main__":
    
    

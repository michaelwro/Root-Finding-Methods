"""

ROOT-FINDING ALGORITHMS: THE SECANT METHOD

Determine the roots of a function via the Secant Method.


Code By:
    Michael Wrona (B.S. Aerospace Engineering)

Find Me At:
    YouTube: @MicWro Engr
    GitHub:  @michaelwro
    Blog:    mwrona.com

"""


import numpy as np


def SecantMethod(funct, p0, p1, TOL=1e-8, MAX_ITER=500, PRINT_RES=False):
    """Determine the root of a function 'funct' between 'a' and 'b' with the Secant Method.

    The function must be continuous and defined on the interval [a, b].

    Args
    ----
        funct (function): Function to compute the root of.
        p0 (float): Left boundary point.
        p1 (float): Right boundary point.
        TOL (float): Solution tolerance. Default 1e-8.
        MAX_ITER (int): Maximum number of allowed iterations. Default 500.
        PRINT_RES (bool): Print solver/solution results. Default False.
    
    Returns
    -------
        soln (float): Computed root.
    """

    p = 0.0  # Store previous solution
    soln = 0.0  # Root solution

    # Check if there's a sign change between the boundary points. If not, the root is not
    # bracketed. f(p0) and f(p1) must have opposite signs
    if (funct(p0) * funct(p1)) >= 0.0:
        print("ERROR IN SECANT_METHOD. NO SIGN CHANGE DETECTED WITH GIVEN BOUNDARY POINTS. SELECT DIFFERENT P0 AND P1.")
        return 0


    for iters in range(MAX_ITER):

        # Approximate the derivative with a secant line
        fp0 = funct(p0)
        fp1 = funct(p1)
        p = p1 - ((p1 - p0) * fp1) / (fp1 - fp0)

        # Check if tolerance is satisfied
        currTol = np.abs((p - p1) / p1)  # Current relative error
        if currTol < TOL:
            soln = p
            break  # Tolerance is satisfied, break
        
        p0 = p1
        p1 = p
    

    # If max. iters was reached, print error message
    if iters >= MAX_ITER - 1:
        print("MAX ITERATIONS REACHED IN SECANT_METHOD. SOLUTION MAY NOT BE VALID.")
        soln = 0.0
    
    # Print solution and solver results
    if PRINT_RES == True:
        print("SECANT METHOD RESULTS:")
        print("    Root: x = %0.12f" % (soln))
        print("    Iterations: %g" % (iters))
        print("    Relative Error: %g\n" % (currTol))
    
    return soln

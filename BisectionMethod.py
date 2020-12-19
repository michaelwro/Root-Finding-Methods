"""

ROOT-FINDING ALGORITHMS: THE BISECTION METHOD

Determine the roots of a function via the Bisection Method.

A full description of the Bisection Method can be found on my
blog site: https://mwrona.com/posts/05-bisection-method/


Code By:
    Michael Wrona (B.S. Aerospace Engineering)

Find Me At:
    YouTube: @MicWro Engr
    GitHub:  @michaelwro
    Blog:    mwrona.com

"""


import numpy as np


def BisectionMethod(funct, a, b, TOL=1e-8, MAX_ITER=500, PRINT_RES=False):
    """Determine a root of a function 'funct' between 'a' and 'b' with the Bisection Method.

    The function must be continuous and defined on the interval [a, b].

    Args
    ----
        funct (function): Function to compute the root of.
        a (float): Left boundary point.
        b (float): Right boundary point.
        TOL (float): Solution tolerance. Default 1e-8.
        MAX_ITER (int): Maximum number of allowed iterations. Default 500.
        PRINT_RES (bool): Print solver/solution results. Default False.
    
    Returns
    -------
        soln (float): Computed root.
    """

    fa = funct(a)  # Evaluate left point
    soln = 0.0  # Root solution

    # Check if there's a sign change between the boundary points. If not, the root is not
    # bracketed. f(a) and f(b) must have opposite signs
    if (funct(a) * funct(b)) >= 0.0:
        print("ERROR IN BISECTION_METHOD. NO SIGN CHANGE DETECTED WITH GIVEN BOUNDARY POINTS. SELECT DIFFERENT A AND B.")
        return 0


    for iters in range(MAX_ITER):
        p = a + (0.5 * (b - a))  # Determine midpoint of the interval
        fp = funct(p)  # Evaluate the midpoint

        # Check if tolerance is satisfied
        currTol = np.abs(0.5 * (b - a))  # Current error
        if fp == 0.0 or currTol < TOL:
            soln = p
            break  # Tolerance is satisfied, break
        
        # Determine new interval to check, i.e. search left or right of the midpoint
        if (np.sign(fa) * np.sign(fp)) > 0.0:
            # If positive, move left
            a = p
            fa = fp
        else:
            # Otherwise, move right
            b = p


    # If max. iters was reached, print error message
    if iters >= MAX_ITER - 1:
        print("MAX ITERATIONS REACHED IN BISECTION_METHOD. SOLUTION MAY NOT BE VALID.")
    
    # Print solution and solver results
    if PRINT_RES == True:
        print("BISECTION METHOD RESULTS:")
        print("    Root: x = %0.12f" % (soln))
        print("    Iterations: %g" % (iters))
        print("    Error: %g\n" % (currTol))
    
    return soln  # Return root solution

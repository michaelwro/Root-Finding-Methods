"""

ROOT-FINDING ALGORITHMS: BRENT'S METHOD

Determine the roots of a function via Brent's Method.


Code By:
    Michael Wrona (B.S. Aerospace Engineering)

Find Me At:
    YouTube: @MicWro Engr
    GitHub:  @michaelwro
    Blog:    mwrona.com

"""


import numpy as np


def BrentsMethod(funct, a, b, TOL=1e-8, MAX_ITER=500, PRINT_RES=False):
    """Determine a root of a function 'funct' between 'a' and 'b' with Brent's Method.

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

    # Check if there's a sign change between the boundary points. If not, the root is not
    # bracketed. f(a) and f(b) must have opposite signs.
    fa = funct(a)
    fb = funct(b)
    if (fa * fb) >= 0.0:
        print("ERROR IN SECANT_METHOD. NO SIGN CHANGE DETECTED WITH GIVEN BOUNDARY POINTS. SELECT DIFFERENT P0 AND P1.")
        return 0
    
    # Swap boundary points if necessary
    if np.abs(fa) < np.abs(fb):
        a, b = b, a
        fa, fb = fb, fa
    
    # Initialize variables
    c = a
    d = a
    soln = 0.0
    bisectFlag = True  # Bisection method flag
    EPS2 = 2.0 * np.abs(np.finfo(float).eps)


    for iters in range(MAX_ITER):

        # Perform tolerance check
        currTol = np.abs(b - a)
        if (currTol < TOL):
            soln = b  # Tolerance met, break loop
            break
        
        # Check if inverse quadratic interpolation needs to be performed
        fc = funct(c)
        
        if np.abs(fa - fc) > EPS2 and np.abs(fb - fc) > EPS2:
            # Inverse quadratic interpolation (Lagrange polynomial)
            term1 = (a * fb * fc) / ((fa - fb) * (fa - fc))
            term2 = (b * fa * fc) / ((fb - fa) * (fb - fc))
            term3 = (c * fa * fb) / ((fc - fa) * (fc - fb))
            p = term1 + term2 + term3
        else:
            p = b - (fb * (b - a) / (fb - fa))  # Otherwise, use the Secant Method
        
        # Check wheter to use the Bisection Method or not
        term1 = 0.25 * ((3.0 * a) + b)
        if (p < term1 and p > b) or \
            (bisectFlag == True and np.abs(p - b) >= (0.5 * np.abs(b - c))) or \
            (bisectFlag == False and np.abs(p - b) >= (0.5 * np.abs(c - d))) or \
            (bisectFlag == True and np.abs(b - c) < EPS2) or \
            (bisectFlag == False and np.abs(c - d) < EPS2):
            p = 0.5 * (a + b)  # Bisection Method
            bisectFlag = True
        else:
            bisectFlag = False  # Otherwise, set bisectFlag to False
        
        fp = funct(p)
        d = c
        c = b

        # If f(a) and f(p) have opposite signs
        if (fa * fp) < 0.0:
            b, fb = p, fp
        else:
            # Otherwise, set a equal to p
            a, fa = p, fp
        
        # If |f(a)| is less than |f(b)|, swap a and b
        if np.abs(fa) < np.abs(fb):
            a, b = b, a
            fa, fb = fb, fa
    

    # If max. iters was reached, print error message
    if iters >= MAX_ITER - 1:
        print("MAX ITERATIONS REACHED IN BRENTS_METHOD. SOLUTION MAY NOT BE VALID.")
    
    # Print solution and solver results
    if PRINT_RES == True:
        print("BRENTS METHOD RESULTS:")
        print("    Root: x = %0.12f" % (soln))
        print("    Iterations: %g" % (iters))
        print("    Error: %g" % (currTol))
    
    return soln  # Return solution

"""

EXAMPLE OF USING THE BISECTION METHOD, SECANT METHOD, AND BRENT'S METHOD
FOR COMPUTING THE ROOT OF A FUNCTION

Code By:
    Michael Wrona (B.S. Aerospace Engineering)

Find Me At:
    YouTube: @MicWro Engr
    GitHub:  @michaelwro
    Blog:    mwrona.com

"""

import numpy as np
import matplotlib.pyplot as plt
from BisectionMethod import BisectionMethod
from SecantMethod import SecantMethod
from BrentsMethod import BrentsMethod



def MyFunct(x):
    """Insert your function here. This is the function we will compute the root of."""
    fx = (-1.0 * (x - 2.0)**2) + 1.0 
    return fx


leftPoint = 1.5  # Left boundary point where the zero is between
rightPoint = 4  # Right boundary point where the zero is between

# SOLVE VIA BISECTION METHOD
bisectionResult = BisectionMethod(MyFunct, leftPoint, rightPoint, PRINT_RES=True)

# SOLVE VIA SECANT METHOD
secantResults = SecantMethod(MyFunct, leftPoint, rightPoint, PRINT_RES=True)

# SOLVE VIA BRENT'S METHOD
brentsResults = BrentsMethod(MyFunct, leftPoint, rightPoint, PRINT_RES=True)



# PLOT THE FUNCTION
x_points = np.linspace(leftPoint, rightPoint, 80)  # Points to pass to function

plt.figure()  # Plot f(x)
plt.plot(x_points, MyFunct(x_points))
plt.title('Function Plot')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.grid()

plt.show()


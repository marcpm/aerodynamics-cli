import argparse
import sympy as sy
import numpy as np
import sys

def tat_cambered(naca, alpha, eqs=None, indefinite=False):
    theta, f, p, t = sy.symbols("theta f p t")
    naca = str(naca)
    f = float(naca[0]) * 0.01
    p = float(naca[1]) * 0.1
    t = float(naca[2]) * 0.01

    if eqs is None:
        lower = 0
        theta_p = np.arccos(1-2*p)
        upper = sy.pi
        zm1 = f/p**2 * (2*p -1 + sy.cos(theta))
        zm2 = f/(1-p)**2 * (2*p -1 + sy.cos(theta))
    else:
        zm1, zm2, lower, theta_p, upper = eqs["zm1"], eqs["zm2"], eqs["lower"], eqs["theta_p"], eqs["upper"]

    A_0 = alpha -  1/sy.pi *( sy.integrate(zm1, (theta, lower,  theta_p)) + sy.integrate(zm2, (theta, theta_p, upper)))
    A_1 = 2/sy.pi * (sy.integrate(zm1* sy.cos(theta), (theta, lower, theta_p)) + sy.integrate(zm2 * sy.cos(theta), (theta, theta_p, upper)))
    A_2 = 2/sy.pi * (sy.integrate(zm1* sy.cos(2*theta), (theta, lower, theta_p)) + sy.integrate(zm2 * sy.cos(2*theta), (theta, theta_p, upper)))
    alpha_l0 = -1/sy.pi * (sy.integrate(zm1 * (sy.cos(theta) -1), (theta, lower, theta_p)) + sy.integrate(zm2 * (sy.cos(theta)-1), (theta, theta_p, upper)))

    cl = sy.pi * (2*A_0 + A_1)
    cm_le = -(A_0 + A_1 - A_2/2)*sy.pi/2
    cm_0 = sy.pi/4 * (A_2 - A_1)
    cm_le2 = -cl/4 + cm_0
    

    print("\nComputations for NACA-{}  at AOA: {} degrees ".format(naca, alpha*180/np.pi))
    if not indefinite:
        print("\nalpha_l0 = {:.5f}".format((180/sy.pi  * alpha_l0).evalf()))
        print("cl = {:.5f}".format(cl.evalf()))
        print("cm_le = {:.5f}".format(cm_le.evalf()))
        print("cm_0 = {:.5f}\n".format(cm_0.evalf()))
    else:
        simpl_func = lambda x: sy.expand(x)
        print("\nalpha_l0 = {}".format(simpl_func(180.0/sy.pi  * alpha_l0)))
        print("cl = {}".format(simpl_func(cl)))
        print("cm_le = {}".format(simpl_func(cm_le)))
        print("cm_0 = {}\n".format(simpl_func(cm_0)))

def main():
    parser = argparse.ArgumentParser(description='Thin Airfoil Theory Computer')
    parser.add_argument('naca',  type=int,
                help='NACA-xxxxxx airfoil')
    parser.add_argument('alpha', type=int,
                help='Angle of Attack', nargs="?", default=0.0)
    parser.add_argument('eq', type=str, help='Mean camber line Equation, by parts allowed', nargs="?", default=None)
    args = parser.parse_args()

    naca = args.naca
    alpha = args.alpha / 180 * np.pi
    eq = args.eq
    # naca, alpha = sys.argv[1] , sys.argv[2]
    tat_cambered(naca, alpha, eq)

def manual():
    theta, eps = sy.symbols("theta  eps")
    naca = 4412     # int NACA
    alpha = 3.0     # degrees
    eqs = None
    eqs = {"zm1": eps * (sy.cos(theta) - 0.5),
            "lower": 0.0, "theta_p": sy.pi/2,
            "zm2": -eps * (sy.cos(theta) + 0.5),
                          "upper": sy.pi}

    indefinite = True

    alpha =  alpha / 180 * np.pi

    # naca, alpha = sys.argv[1] , sys.argv[2]
    tat_cambered(naca, alpha, eqs, indefinite)
if __name__ == "__main__":

    manual()


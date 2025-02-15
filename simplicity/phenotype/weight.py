import numpy as np
from scipy.optimize import fsolve

t_max = 21

def w_t_params(t_half=30, t_max=t_max):
    """
    Compute values for k_e and k_a from known pharmacokinetics relationships

    Parameters
    ----------
    t_half : int, optional
        half life of antibodies in plasma. The default is 30.
    t_max : int, optional
        time of max antibody concentration after exposure. The default is t_max.

    Returns
    -------
    float
        k_e antibody elimination rate constant.
    float
        k_a - antibody generation rate constant.

    """
    k_e = np.log(2)/t_half # elimination rate constant
    
    # use fsolve to find ka numerically
    
    # Define k_a function to solve it numerically
    def get_k_a(k_a, t_max, k_e):
        # This is the equation derived from t_max = ln(ka/ke) / (ka - ke)
        # We rearrange it to the form f(ka) = 0
        return t_max * (k_a - k_e) - np.log(k_a / k_e)
    
    initial_guess = 0.2  # This is an initial guess for k_a
    # Find the root of the equation
    k_a_solution = fsolve(get_k_a, initial_guess, args=(t_max, k_e))
    k_a = k_a_solution[0] # antibody generation rate constant
    return (k_e,k_a)

# weight function derived from SARS-CoV-2 antibodies pharmacokinetics 
def weights(t,t_eval, k_e, k_a, t_max):
    return ((
            np.exp(k_e * (t-t_eval)) - np.exp(k_a * ((t-t_eval)))) / (
            (np.exp(-k_e * t_max) - np.exp(-k_a * t_max))))


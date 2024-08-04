"""
Calculate Hill Equation coefficients from dose and response values.

Special thanks to Giuseppe Cardillo - https://it.mathworks.com/matlabcentral/fileexchange/38122

This will calculate coefficients with a 95% confidence interval. This is used for dose response curves, ELISA, and other biological techniques.
These are the paremeters used in a Four Parameter Logistic (4PL) Curve, and should match the values of:
"Quest Graph™ Four Parameter Logistic (4PL) Curve Calculator." AAT Bioquest, Inc., 4 Aug. 2024, https://www.aatbio.com/tools/four-parameter-logistic-4pl-curve-regression-online-calculator.

Definitions:


A (Minimum Asymptote): The minimum response value (or asymptote) that the curve approaches as the independent variable (e.g., concentration) approaches negative infinity.
In a bioassay or dose-response curve, this can be thought of as the response value at zero standard concentration. Essentially, it’s the lower plateau of the sigmoidal curve.

B (Hill’s Slope): Refers to the steepness of the curve. Positive values of B make the curve rise (sigmoidal shape), while negative values make it fall.
The Hill slope characterizes how quickly the response changes with increasing concentration.
A higher absolute value of B indicates a steeper curve.

C (Inflection Point):
Represents the concentration of analyte (e.g., ligand, protein) where the curve changes direction or signs.
At C, the response is halfway between the minimum and maximum asymptotes (i.e., y = (D - A) / 2).
The inflection point is where the curvature transitions from concave down to concave up (or vice versa).

D (Maximum Asymptote): Represents the maximum 
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def L4P(x, A, B, C, D):
    """
    Four Parameters Logistic Regression (4PL) model.
    F(x) = D + (A - D) / (1 + (x / C)**B)
    """
    return D + (A - D) / (1 + (x / C)**B)


def calculate_values(concentration_array, response_array, curve_direction):

    # Flip order of both arrays, if concentration is now smallest to highest
    if concentration_array[0]>concentration_array[-1]:
        concentration_array.reverse()
        response_array.reverse()

    if concentration_array[0] == 0:
        # Log of zero is invalid, so set it to a near zero value
        concentration_array[0] = 0.000000000000000000000001
    x_data = np.array(concentration_array)
    #x_data = np.log10(x_data)
    y_data = np.array(response_array)

    # Initial parameter guesses (you can adjust these)
    if curve_direction == "up":
        # For curves that go up
        guess_neg_inf_asymptote = 0.001
        guess_hill_slope = 1.515
        guess_ec50 = 108
        guess_pos_inf_asymptote = 3.784
    else:
        # For curves that go down
        guess_neg_inf_asymptote = 10
        guess_hill_slope = -0.3
        guess_ec50 = 0.4
        guess_pos_inf_asymptote = 90

    initial_guess = [guess_neg_inf_asymptote, guess_hill_slope, guess_ec50, guess_pos_inf_asymptote]

    # Fit the 4PL model to the data
    params, covariance = curve_fit(L4P, x_data, y_data, p0=initial_guess, maxfev=3000000)

    # Extract the fitted parameters
    A_fit, B_fit, C_fit, D_fit = params

    if np.isnan(A_fit):
        print(f"Estimated parameters: A = {A_fit:.6f}, B = {B_fit:.6f}, C = {C_fit:.6f}, D = {D_fit:.6f}")
        raise Exception("Failure to derive values")

    # Calculate goodness-of-fit measures (you can customize this)
    residuals = y_data - L4P(x_data, *params)
    sse = np.sum(residuals**2)
    rmse = np.sqrt(sse / len(x_data))
    r2 = 1 - sse / np.sum((y_data - np.mean(y_data))**2)

    return A_fit, B_fit, C_fit, D_fit, params, x_data, y_data, sse, rmse, r2


def print_results(A_fit, B_fit, C_fit, D_fit, sse, rmse, r2, curve_direction, is_final=False):
    print("\n\n######################################")
    # Print the results
    if not is_final:
        print("Using a test curve direction of: "+curve_direction)
        print(f"Estimated parameters: A = {A_fit:.6f}, B = {B_fit:.6f}, C = {C_fit:.6f}, D = {D_fit:.6f}")
        print(f"Goodness-of-fit measures: SSE = {sse:.6f}, RMSE = {rmse:.6f}, R^2 = {r2:.6f}")
    else:
        print("Final curve direction: "+curve_direction)
        print(f"Final parameters: A = {A_fit:.6f}, B = {B_fit:.6f}, C = {C_fit:.6f}, D = {D_fit:.6f}")
        print(f"Goodness-of-fit measures: SSE = {sse:.6f}, RMSE = {rmse:.6f}, R^2 = {r2:.6f}")

def render_plot(x_data, y_data, params):
    plt.scatter(x_data, y_data, label="Data", color="red")
    x_fit = np.linspace(min(x_data), max(x_data), 100)
    y_fit = L4P(x_fit, *params)
    plt.plot(x_fit, y_fit, label="Fitted Curve", color="blue")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

def calculate_coefficients(concentration_array, response_array):
    A_fit_final = None
    B_fit_final = None
    C_fit_final = None
    D_fit_final = None
    params_final = None
    x_data_final = None
    y_data_final = None
    sse_final = None
    rmse_final = None
    r2_final = None
    curve_direction_final = None
    
    # Try one curve direction first, but store the best R^2 as the final
    for possible_curve_direction in ["up","down"]:
        A_fit, B_fit, C_fit, D_fit, params, x_data, y_data, sse, rmse, r2 = calculate_values(concentration_array, response_array, possible_curve_direction)
        print_results(A_fit, B_fit, C_fit, D_fit, sse, rmse, r2, possible_curve_direction)

        # Only pick the direction of curve that has the highest R-squared value
        if r2_final is None or r2 > r2_final:
            A_fit_final = A_fit
            B_fit_final = B_fit
            C_fit_final = C_fit
            D_fit_final = D_fit
            params_final = params
            x_data_final = x_data
            y_data_final = y_data
            sse_final = sse
            rmse_final = rmse
            r2_final = r2
            curve_direction_final = possible_curve_direction
    return A_fit_final, B_fit_final, C_fit_final, D_fit_final, params_final, x_data_final, y_data_final, sse_final, rmse_final, r2_final, curve_direction_final


if __name__ == "__main__":

    """
    Example Data #1 -

    Source: 
    https://www.medcalc.org/manual/nonlinearregression-example.php
    # See also this useful calculator - https://www.aatbio.com/tools/four-parameter-logistic-4pl-curve-regression-online-calculator

    Expected output:
    A = 0.1536
    B = 1.7718
    C = 19.3494
    D = 28.4479
    """
    # Example 1 data (replace with your actual data)
    concentration_array = [0, 1.3, 2.8, 5, 10.2, 16.5, 21.3, 31.8, 52.2]
    response_array = [0.1, 0.5, 0.9, 2.6, 7.1, 12.3, 15.3, 20.4, 24.2]

    """
    Example Data #2
    
    Source: 
    Cardillo G. (2012) Four parameters logistic regression - There and back again
    https://it.mathworks.com/matlabcentral/fileexchange/38122

    Coefficients (with 95% confidence bounds):
    A =    0.001002  (-0.04594, 0.04794)
    B =       1.515  (1.293, 1.738)
    C =         108  (86.58, 129.4)
    D =       3.784  (3.302, 4.266)

    sse: 0.0012
    rsquare: 0.9998
    dfe: 3
    adjrsquare: 0.9996
    rmse: 0.0200
    """
    # Example 2 data
    # concentration_array = [0, 4.5, 10.6, 19.7, 40, 84, 210]
    # response_array = [0.0089, 0.0419, 0.0873, 0.2599, 0.7074, 1.528, 2.7739]

    A_fit_final, B_fit_final, C_fit_final, D_fit_final, params_final, x_data_final, y_data_final, sse_final, rmse_final, r2_final, curve_direction_final = calculate_coefficients(concentration_array, response_array)

    print("Using "+curve_direction_final+" version of curve data")
    print_results(A_fit_final, B_fit_final, C_fit_final, D_fit_final, sse_final, rmse_final, r2_final, curve_direction_final, is_final=True)
    render_plot(x_data_final, y_data_final, params_final)

   
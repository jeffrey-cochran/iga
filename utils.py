from numpy import searchsorted, s_, broadcast_to, vstack, ones, zeros, atleast_2d
from numpy.typing import ArrayLike, NDArray

def find_span(parameter_values: NDArray, knot_vector: NDArray) -> NDArray:
    """Algorithm A2.1 from Piegl and Tiller, extended to a vector of paramter values"""
    if parameter_values.min() < knot_vector.min():
        raise ValueError("Parameter value cannot be less than the least knot.")
    elif parameter_values.max() > knot_vector.max():
        raise ValueError("Parameter value cannot be greater than the greatest knot.")

    precedes = searchsorted(knot_vector, parameter_values, side='left')
    succeeds = precedes - 1
    succeeds[succeeds < 0] = 0

    return succeeds


def basis_funcs(knot_spans: NDArray, parameter_values: NDArray, knot_vector: NDArray, poly_order: int, verbose: bool = False):
    """Algorithm A2.2 from Piegl and Tiller, extended to a vector of parameter values"""

    if knot_spans.shape != parameter_values.shape:
        raise ValueError("The number of parameter values must be the same as the number of corresponding knot spans, but: knot_spans.shape != parameter_values.shape")

    #
    # Compute the slices 
    left_slices = tuple((
        s_[i:j:-1] for (i,j) in zip(knot_spans, knot_spans-poly_order )
    )) 
    right_slices = tuple((
        s_[i:j:1] for (i,j) in zip(knot_spans+1, knot_spans+poly_order+1)
    ))

    #
    # Compute knot differences 
    left = parameter_values - vstack( 
        tuple((
            knot_vector[l_] for l_ in left_slices
        ))
    ).transpose()
    right = vstack(
        tuple((
            knot_vector[r_] for r_ in right_slices
        ))
    ).transpose() - parameter_values

    #
    # Instantiate output
    out_shape = (poly_order+1, poly_order+1, parameter_values.shape[0]) if verbose else (poly_order+1, parameter_values.shape[0]) 
    func_values = ones(out_shape) 

    #
    # Iterate over polynomial orders
    for j in range(1, poly_order+1):
        saved = zeros((parameter_values.shape[0],))
        _slice = s_[j,j,:] if verbose else s_[j,:]
        for r in range(j):
            denom = right[r,:] + left[j-r-1,:]
            if verbose:
                #
                # lower triangle
                func_values[j,r,:] = denom
                temp = func_values[r,j-1,:]/func_values[j,r,:]
                #
                # upper triangle
                func_values[r,j,:] = saved + right[r,:]*temp
                saved = left[j-r-1,:]*temp
            else:
                temp = func_values[r,:] / denom
                func_values[r,:] = saved + right[r,:]*temp
                saved = left[j-r-1, :]*temp
        func_values[_slice] = saved
    
    return func_values

def ders_basis_funcs(knot_spans: NDArray, parameter_values: NDArray, knot_vector: NDArray, poly_order: int, der_order: int) -> NDArray:
    """Algorithm A2.3 from Piegl and Tiller, extended to a vector of parameter values"""

    #
    # Instantiate output
    num_paramater_values = parameter_values.shape[0]
    ders = zeros((der_order + 1, poly_order + 1, num_paramater_values))
    #
    # Compute basis function values and knot differences
    ndu = basis_funcs(knot_spans, parameter_values, knot_vector, poly_order, verbose=True)

    #
    # Set 0th derivative to value of the function (pth order polynomial)
    ders[0,:,:] = ndu[:,poly_order,:]

    #
    # Loop over function index
    for func_idx in range(poly_order+1):
        s1 = 0
        s2 = 1
        temp = zeros((2, der_order+1, num_paramater_values))
        temp[0,0,:] = 1
        for der_idx in range(1, der_order+1):
            d = zeros((num_paramater_values,))
            rk = func_idx - der_idx
            pk = poly_order - der_idx

            if rk >= 0:
                temp[s2,0,:] = temp[s1,0,:] / ndu[pk+1,rk,:]
                d = temp[s2,0,:]*ndu[rk,pk,:]
            
            j1 = 1 if rk >= -1 else -rk
            j2 = der_idx-1 if (func_idx-1) <= pk else (poly_order - func_idx)
            for j in range(j1, j2+1):
                temp[s2,j,:] = (temp[s1,j,:] - temp[s1, j-1,:]) / ndu[pk+1,rk+j,:]
                d += temp[s2,j,:]*ndu[rk+j,pk,:]

            if func_idx <= pk:
                temp[s2,der_idx,:] = -temp[s1,der_idx-1,:] / ndu[pk+1,func_idx,:]
                d += temp[s2,der_idx,:]*ndu[func_idx,pk,:]

            ders[der_idx,func_idx,:] = d
            
            #
            # Swap indices
            dummy = s1
            s1 = s2
            s2 = dummy

        mult_factor = poly_order
        for der_idx in range(1,der_order+1):
            ders[der_idx,:,:] *= mult_factor
            mult_factor *= (poly_order-der_idx)

    return ders

def curve_point(parameter_values: NDArray, knot_vector: NDArray, control_points: NDArray, poly_order) -> NDArray:
    span_eval = find_span(parameter_values, knot_vector)
    basis_func_eval = basis_funcs(span_eval, parameter_values, knot_vector, poly_order, verbose=False)

    curve_points = zeros((parameter_values.shape[0], control_points.shape[1]))
    for i in range(poly_order+1):
        curve_points += (basis_func_eval[i,:]*(control_points[span_eval-poly_order+i, :]).transpose()).transpose()
        print(curve_points)

    return curve_points

def plot_bspline(knot_vector: NDArray, control_points: NDArray, poly_order: int) -> NDArray:
    unique_knots = 
    return
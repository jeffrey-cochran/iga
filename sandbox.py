from constants import test_knots, test_parameter_values
from utils import find_span, basis_funcs

knot_span = find_span(test_parameter_values, test_knots)

print(
    basis_funcs(knot_span, test_parameter_values, test_knots, 2, verbose=True)
)
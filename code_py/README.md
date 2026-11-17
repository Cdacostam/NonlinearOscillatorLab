# nonlinear_oscillator_lab (Python package, v0.2.2)

Energy geometry, periods (double and arbitrary precision), closed-form solutions and simulation of
x'' + a x + b x^3 + c x^5 = 0, plus a damped, harmonically forced extension.

From the bundle root:

    python -m pip install './code_py[test]'
    python -m pytest code_py/tests -q
    python code_py/examples/quickstart.py

Documentation: docs/ (installation, tutorial, API, troubleshooting, reproducibility, scope).
The regression tests read ../data/wolfram/; the benchmark test runs ../reproduce/python_benchmarks.py.

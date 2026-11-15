# nonlinear_oscillator_lab (Python package, v0.2.1)

Energy geometry, periods (double and arbitrary precision), closed-form solutions and simulation of
x'' + a x + b x^3 + c x^5 = 0, plus a damped, harmonically forced extension.

    python -m pip install .            # from this folder
    python -m pytest -q                # 166 tests: pytest -m unit|regression|failure|integration
    python examples/quickstart.py

Documentation: docs/ (installation, tutorial, API, troubleshooting, reproducibility, scope).
The regression tests read ../data/wolfram/; the benchmark test runs ../reproduce/python_benchmarks.py.

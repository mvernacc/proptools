proptools: Rocket Propulsion Design Tools
=========================================

Proptools is a Python package for preliminary design of rocket propulsion systems. 

Proptools provides implementations of equations for nozzle flow, turbo-machinery and rocket structures.
The project aims to cover most of the commonly used equations in
*Rocket Propulsion Elements* by George Sutton and Oscar Biblarz and
*Modern Engineering for Design of Liquid-Propellant Rocket Engines* by Dieter Huzel and David Huang.

Proptools can be used as a desktop calculator:

.. code-block:: python

    >> from proptools import nozzle
    >> p_c = 10e6; p_e = 100e3; gamma = 1.2; m_molar = 20e-3; T_c = 3000.
    >> C_f = nozzle.thrust_coef(p_c, p_e, gamma)
    >> c_star = nozzle.c_star(gamma, m_molar, T_c)
    >> I_sp = C_f * c_star / nozzle.g
    >> print "The engine's ideal sea level specific impulse is {:.1f} seconds.".format(I_sp)
    The engine's ideal sea level specific impulse is 288.7 seconds.

Proptools can also be used as a library in other propulsion design and analysis software. It is distributed under a `MIT License <https://github.com/mvernacc/proptools/blob/master/LICENSE>`_ and can be used in commercial projects.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   modules



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

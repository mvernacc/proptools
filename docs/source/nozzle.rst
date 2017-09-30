Nozzle Flow
***********

The purpose of a rocket is to generate thrust by expelling mass at high velocity. A nozzle is a flow device which accelerates gas to high velocity before it is expelled from the vehicle. The nozzle accelerates the gas by converting some of the gas's thermal energy into kinetic energy.


Ideal Nozzle Flow
=================
Ideal nozzle flow is a simplified model of the aero- and thermo-dynamic behavior of fluid in a nozzle. The ideal model allows us to write algebraic relations between an engine's geometry and operating conditions (e.g. throat area, chamber pressure, chamber temperature) and its performance (e.g. thrust and specific impulse). These equations are fundamental tools for the preliminary design of rocket propulsion systems.


The assumptions of the ideal model are:

#. The fluid flowing through the nozzle is gaseous.
#. The gas is homogeneous, obeys the ideal gas law, and is calorically perfect. Its molar mass (:math:`\mathcal{M}`) and heat capacities (:math:`c_p, c_v`) are constant throughout the fluid, and do not vary with temperature.
#. There is no heat transfer to or from the gas. Therefore, the flow is adiabatic. The specific enthalpy :math:`h` is constant throughout the nozzle.
#. There are no viscous effects or shocks within the gas or at its boundaries. Therefore, the flow is reversible. If the flow is both adiabatic and reversible, it is isentropic: the specific entropy :math:`s` is constant throughout the nozzle.
#. The flow is steady; :math:`\frac{d}{dt} = 0`.
#. The flow is quasi one dimensional. The flow state varies only in the axial direction of the nozzle.
#. The flow velocity is axially directed.
#. The flow does not react in the nozzle. The chemical equilibrium established in the combustion chamber does not change as the gas flows through the nozzle. This assumption is known as "frozen flow".

These assumptions are usually acceptably accurate for preliminary design work. Most rocket engines perform within 1% to 6% of the ideal model predictions [RPE]_.


Isentropic Relations
====================
Under the assumption of isentropic flow and calorically perfect gas, there are several useful relations between fluid states. These relations depend on the heat capacity ratio, :math:`\gamma = c_p /c_v`. Consider two gas states, 1 and 2, which are isentropically related (:math:`s_1 = s_2`). The states' pressure, temperature and density ratios are related:

.. math::

  \frac{p_1}{p_2} = \left( \frac{\rho_1}{\rho_2} \right)^\gamma = \left( \frac{T_1}{T_2} \right)^{\frac{\gamma}{\gamma - 1}}

Stagnation state
----------------

Now consider the relation between static and stagnation states in a moving fluid. The stagnation state is the state a moving fluid would reach if it were isentropically decelerated to zero velocity. The stagnation enthalpy :math:`h_0` is the sum of the static enthalpy and the specific kinetic energy:

.. math::

  h_0 = h + \frac{1}{2} v^2

For a calorically perfect gas, :math:`T = h / c_p`, and the stagnation temperature is:

.. math::

  T_0 = T + \frac{v^2}{2 c_p}

It is helpful to write the fluid properties in terms of the Mach number :math:`M`, instead of the velocity. Mach number is the velocity normalized by the local speed of sound, :math:`a = \sqrt{\gamma R T}`. In terms of Mach number, the stagnation temperature is:

.. math::
  T_0 = T \left( 1 + \frac{\gamma - 1}{2} M^2 \right)

Because the static and stagnation states are isentropically related, :math:`\frac{p_0}{p} = \left( \frac{T_0}{T} \right)^{\frac{\gamma}{\gamma - 1}}`. Therefore, the stagnation pressure is:

.. math::
  p_0 = p \left( 1 + \frac{\gamma - 1}{2} M^2 \right)^{\frac{\gamma}{\gamma - 1}}

Use ``proptools`` to plot the stagnation state variables against Mach number:

TODO code example

Exit velocity
-------------

The exit velocity of the exhaust gas is the fundamental measure of efficiency for rocket propulsion systems, as the `rocket equation <https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation>`_ shows. We can now show a basic relation between the exit velocity and the combustion conditions of the rocket. First, use the conservation of energy to relate the velocity at any two points in the flow:

.. math::

  v_2 = \sqrt{2(h_1 - h_2) + v_1^2}

We can replace the enthalpy difference with an expression of the pressures and temperatures, using the isentropic relations.

.. math::

  v_2 = \sqrt{\frac{2 \gamma}{\gamma - 1} R T_1 \left(1 - \left( \frac{p_2}{p_1} \right)^{\frac{\gamma - 1}{\gamma}} \right) + v_1^2}

Set state 1 to be the conditions in the combustion chamber: :math:`T_1 = T_c, p_1 = p_c, v_1 \approx 0`. Set state 2 to be the state at the nozzle exit: :math:`p_2 = p_e, v_2 = v_e`. This gives the exit velocity:

.. math::

  v_e &= \sqrt{\frac{2 \gamma}{\gamma - 1} R T_c \left(1 - \left( \frac{p_e}{p_c} \right)^{\frac{\gamma - 1}{\gamma}} \right)} \\
  &= \sqrt{\frac{2 \gamma}{\gamma - 1} \mathcal{R} \frac{T_c}{\mathcal{M}} \left(1 - \left( \frac{p_e}{p_c} \right)^{\frac{\gamma - 1}{\gamma}} \right)}

.. math::

  \require{siunitx}

where :math:`\mathcal{R} = 8.314` J mol :sup:`-1` K :sup:`-1` is the universal gas constant and :math:`\mathcal{M}` is the molar mass of the exhaust gas.



TODO velocity from energy conservation


Mach-area Relation
==================
TODO need for converging-diverging nozzle
TOOD importance of expansion ratio

Choked Flow
===========
TODO

Thrust
======
TODO

Thrust coefficient
==================
TODO

Characteristic velocity
=======================
TODO

Specific Impulse
================
TODO



.. [RPE] G. P. Sutton and O. Biblarz, *Rocket Propulsion Elements*, Hoboken: John Wiley & Sons, 2010.
